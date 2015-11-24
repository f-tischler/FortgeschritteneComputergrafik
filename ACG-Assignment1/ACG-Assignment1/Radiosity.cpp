/******************************************************************
*
* Radiosity.cpp
*
* Description: This file demonstrates global illumination rendering
* based on the radiosity method. The geometry is divided into patches
* for which the form factors are determined employing Monte Carlo 
* integration. Radiosity values for the patches are computed with
* an iterative solver. The final image (i.e. radiance) is obtained 
* via tracing rays into the scene. Two output files are saved - 
* one with constant shading of patches and one with bicubic color
* interpolation.
*
* The code is extended from software by user Hole and Kevin Beason,
* released under the MIT License. 
*
* http://kagamin.net/hole/license.txt
* http://kagamin.net/hole/smallpt-license.txt
*
* Advanced Computer Graphics Proseminar WS 2015
* 
* Interactive Graphics and Simulation Group
* Institute of Computer Science
* University of Innsbruck
*
*******************************************************************/

#ifdef _WIN32

#define _CRT_SECURE_NO_WARNINGS

#include <math.h>

#define _USE_MATH_DEFINES

#define drand48() (((double)rand())/((double)RAND_MAX))

#endif

/* Standard includes */
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include <assert.h>

using namespace std;

const double Over_M_PI = 1.0/M_PI;

static double *form_factor;
static int patch_num = 0;

/*------------------------------------------------------------------
| Struct for standard vector operations in 3D 
| (used for points, vectors, and colors)
------------------------------------------------------------------*/
struct Vector 
{
    double x, y, z;           /* Position XYZ or color RGB */

    Vector(const Vector &b) : x(b.x), y(b.y), z(b.z) {}
    Vector(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    
    Vector operator+(const Vector &b) const 
    {
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator-(const Vector &b) const
    {
        return Vector(x - b.x, y - b.y, z - b.z);
    }

    Vector operator/(double c) const 
    {
        return Vector(x / c, y / c, z / c);
    }

    Vector operator*(double c) const 
    {
        return Vector(x * c, y * c, z * c);
    }

    friend Vector operator*(double c, const Vector &b) 
    { 
        return b * c; 
    }

    Vector MultComponents(const Vector &b) const
    {
        return Vector(x * b.x, y * b.y, z * b.z);
    }

    const double LengthSquared() const 
    { 
        return x*x + y*y + z*z; 
    }

    const double Length() const 
    { 
        return sqrt(LengthSquared()); 
    } 

    const Vector Normalized() const
    {
        return Vector(x, y, z) / sqrt(x*x + y*y + z*z);
    }

    const double Dot(const Vector &b) const 
    {
        return x * b.x + y * b.y + z * b.z;
    }

    const Vector Cross(const Vector &b) const
    {
        return Vector((y * b.z) - (z * b.y), 
                      (z * b.x) - (x * b.z), 
                      (x * b.y) - (y * b.x));
    }
};

typedef Vector Color;
const Color BackgroundColor(0.0, 0.0, 0.0);

/*------------------------------------------------------------------
| Structure for rays (e.g. viewing ray, ray tracing)
------------------------------------------------------------------*/
struct Ray 
{
	Vector org, dir;    /* Origin and direction */
	Ray(const Vector org_, const Vector &dir_) : org(org_), dir(dir_) {}
};

/*------------------------------------------------------------------
| Struct holds pixels/colors of rendered image
------------------------------------------------------------------*/
struct Image 
{
    int width, height;
    Color *pixels;

    Image(int _w, int _h) : width(_w), height(_h) 
    {
        pixels = new Color[width * height];        
    }

    Color getColor(int x, int y) 
    {
        int image_index = (height-y-1) * width + x;	
        return pixels[image_index];
    }

    void setColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x;	
        pixels[image_index] = c;
    }

    void addColor(int x, int y, const Color &c) 
    {
        int image_index = (height-y-1) * width + x;	
        pixels[image_index] = pixels[image_index] + c;
    }

    int toInteger(double x)
    { 
        /* Clamp to [0,1] */
        if (x<0.0)
            x = 0.0;		
        
        if (x>1.0)
            x = 1.0;             

        /* Apply gamma correction and convert to integer */
        return int(pow(x,1/2.2)*255+.5); 
    }

    void Save(const string &filename) 
    {
        /* Save image in PPM format */
        FILE *f = fopen(filename.c_str(), "wb");
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for (int i = 0; i < width * height; i++)
            fprintf(f,"%d %d %d ", toInteger(pixels[i].x), 
                                   toInteger(pixels[i].y), 
                                   toInteger(pixels[i].z));
        fclose(f);
    }
};

/*------------------------------------------------------------------
| Basic geometric element of scene description;
| Rectangles are subdivided into smaller patches for radiosity
| computation (subdivision equal for all rectangles)
------------------------------------------------------------------*/
struct Rectangle 
{
    Vector p0;                     
    Vector edge_a, edge_b;
    Color emission, color;
    Vector normal;
	
    vector<Color> patch; 
    int a_num, b_num;       /* Number of patches/subdivision of edges */
    double a_len, b_len;    /* Edge lengths */

    Rectangle(const Vector p0_, const Vector &a_, const Vector &b_, 
              const Color &emission_, const Color &color_) :
              p0(p0_), edge_a(a_), edge_b(b_), emission(emission_), color(color_) 
    {
        normal = edge_a.Cross(edge_b);
        normal = normal.Normalized();        
        a_len = edge_a.Length();
        b_len = edge_b.Length();
    }

    Color sample_patch(int ia, int ib) const 
    {
        if (ia < 0) ia = 0;
        if (ia >= a_num) ia = a_num - 1;
        if (ib < 0) ib = 0;
        if (ib >= b_num) ib = b_num - 1;
        return patch[ia * b_num + ib];
    }

    void init_patchs(const int a_num_, const int b_num_) 
    {
        a_num = a_num_;
        b_num = b_num_;
        patch.clear();
        patch.resize(a_num * b_num);
    }

    /* Rectangle-ray intersection */
    const double intersect(const Ray &ray) 
    {
        /* Check for plane-ray intersection first */
        const double t = (p0 - ray.org).Dot(normal) / ray.dir.Dot(normal);
        if (t <= 0.00001)
            return 0.0;

        /* Determine if intersection is within rectangle */
        Vector p = ray.org + ray.dir * t;
        Vector d = p - p0;
        const double ddota = d.Dot(edge_a);
        if (ddota < 0.0 || ddota > edge_a.LengthSquared())
            return 0.0;
        
        const double ddotb = d.Dot(edge_b);
        if (ddotb < 0.0 || ddotb > edge_b.LengthSquared())
            return 0.0;
        
        return t;
    }  
};

class Triangle
{
public:
	Triangle() : _a(0.0) { }

    Triangle(const Vector p0_, const Vector p1_, const Vector p2_, const Color &emission_, const Color &color_) 
		: _p0(p0_), _p1(p1_), _p2(p2_), _emission(emission_), _color(color_)
	{
        _normal = _p0.Cross(_p1).Normalized();
		_a = (_p1 - _p0).Cross(_p2 - _p0).Length();
    }

	double area() const { return _a; }
	Vector normal() const { return _normal; }


	bool intersects(const Ray &ray, float& distance) const 
	{
		Vector e1, e2;  //Edge1, Edge2
		Vector P, Q, T;
		float det, inv_det, u, v;
		float t;

		auto V1 = _p0;
		auto V2 = _p1;
		auto V3 = _p2;
		auto D = ray.dir;
		auto O = ray.org;

		#define SUB(a,b,c) a = (b) - (c)
		#define CROSS(a,b,c) a = (b).Cross(c)
		#define DOT(a,b) (a).Dot(b)
		#define EPSILON 0.000001f

		//Find vectors for two edges sharing V1
		e1 = V2 - V1;
		e2 = V3 - V1;

		//Begin calculating determinant - also used to calculate u parameter
		P = D.Cross(e2);

		//if determinant is near zero, ray lies in plane of triangle
		det = e1.Dot(P);

		//NOT CULLING
		if (det > -EPSILON && det < EPSILON) 
			return false;


		inv_det = 1.f / det;
		//calculate distance from V1 to ray origin
		T = O - V1;

		//Calculate u parameter and test bound
		u = T.Dot(P) * inv_det;

		//The intersection lies outside of the triangle
		if (u < 0.f || u > 1.f) 
			return false;

		//Prepare to test v parameter
		Q = T.Cross(e1);

		//Calculate V parameter and test bound
		v = D.Dot(Q) * inv_det;

		//The intersection lies outside of the triangle
		if (v < 0.f || u + v  > 1.f) 
			return false;

		t = e2.Dot(Q) * inv_det;

		//ray intersection
		if (t > EPSILON) 
		{ 
			distance = t;
			return true;
		}

		return false;
	}

	void split(Triangle(&out)[4]) const
	{
		auto midP1P0 = (_p1 - _p0) / 2;
		auto midP2P0 = (_p2 - _p0) / 2;
		auto midP2P1 = (_p2 - _p1) / 2;

		out[0] = Triangle(midP1P0, midP2P1, midP2P0, _emission, _color);
		out[1] = Triangle(_p0, midP1P0, midP2P1, _emission, _color);
		out[2] = Triangle(midP1P0, _p1, midP2P1, _emission, _color);
		out[3] = Triangle(midP2P1, _p2, midP2P0, _emission, _color);
	}

	void split_n(size_t n, Triangle* out) const
    {
		assert(out);

		for (auto i = 0u; i < n; i++)
		{
			Triangle tmp[4];

			split(tmp);

			for (auto j = 0u; j < 4u; j++)
			{
				out[i * 4 + j] = tmp[j];
			}
		}
    }

	Vector point_inside() const
	{
		auto midP1P0 = (_p1 - _p0) / 2;
		auto midP2P1 = (_p2 - _p1) / 2;

		return (midP2P1 - midP1P0) / 2;
	}

	const Color& getEmission() const { return _emission; }
	const Color& getColor() const { return _color; }
	const Vector& getNormal() const { return _normal; }
	const Vector& getP1() const { return _p0; }
	const Vector& getP2() const { return _p1; }
	const Vector& getP3() const { return _p2; }

private:
	Vector _p0, _p1, _p2, _normal;
	Color _emission, _color;
	double _a;

};

/******************************************************************
* Hard-coded scene definition: the geometry is composed of rectangles. 
* These are defined by:
* - vector to corner(origin), edge a, edge b 
* - emitted light energy (light sources), surface reflectivity (~color)
*******************************************************************/
Rectangle recs[] = 
{	
    /* Cornell Box walls */
    Rectangle(Vector(0.0, 0.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 80.0, 0.0),   
              Vector(), Color(0.75, 0, 0)), /* Back */
    Rectangle(Vector(0.0, 0.0, 170.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, -170.0), 
              Vector(), Color(0, 0.75, 0)), /* Bottom */
    Rectangle(Vector(0.0, 80.0, 0.0), Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0),  
              Vector(), Color(0, 0.75, 0)), /* Top */
    Rectangle(Vector(0.0, 0.0, 170.0), Vector(0.0, 0.0, -170.0), Vector(0.0, 80.0, 0.0),   
              Vector(), Color(0.75, 0.25, 0.25)), /* Left */
    Rectangle(Vector(100.0, 0.0, 0.0), Vector(0.0, 0.0, 170.0), Vector(0.0, 80.0, 0.0),   
              Vector(), Color(0.25, 0.25, 0.75)), /* Right */
    Rectangle(Vector(100.0, 0.0, 170.0), Vector(-100.0, 0.0, 0.0), Vector(0.0, -80.0, 0.0),  
              Vector(), Color(0,1,0)), /* Front (not visible) */
	
    /* Area light source on top */
    Rectangle(Vector(40.0, 79.99, 65.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, 20.0), 
              Vector(12,12,12), Color(0.75, 0.75, 0.75)), 

    /* Cuboid in room */
    Rectangle(Vector(30.0, 0.0, 100.0), Vector(0.0, 0.0, -20.0), Vector(0.0, 40.0, 0.0),   
		Vector(), Color(0.1, 0.1, 0.1)), /* Right */
    Rectangle(Vector(10.0, 0.0, 80.0), Vector(0.0, 0.0, 20.0), Vector(0.0, 40.0, 0.0),   
		Vector(), Color(0.1, 0.1, 0.1)), /* Left */
    Rectangle(Vector(10.0, 0.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 40.0, 0.0),   
		Vector(), Color(0.1, 0.1, 0.1)), /* Front */
    Rectangle(Vector(30.0, 0.0, 80.0), Vector(-20.0, 0.0, 0.0), Vector(0.0, -40.0, 0.0),  
		Vector(), Color(0.1, 0.1, 0.1)), /* Back */
    Rectangle(Vector(10.0, 40.0, 100.0), Vector(20.0, 0.0, 0.0), Vector(0.0, 0.0, -20.0), 
		Vector(), Color(0.1, 0.1, 0.1)), /* Top */
};

std::vector<Triangle> triangles;

template<size_t N>
void convertRectangles(const Rectangle(&rects)[N], std::vector<Triangle>& triangles)
{
	if (!triangles.empty())
		triangles.clear();

	int i = 0;
	triangles.resize(N*2);
	for (auto rect : rects)
	{
		auto tl = rect.p0;
		auto bl = rect.p0 + rect.edge_b;
		auto tr = rect.p0 + rect.edge_a;
		auto br = bl + rect.edge_a;

		triangles[i++] = Triangle(tl, bl, tr, rect.emission, rect.color);
		triangles[i++] = Triangle(bl, tr, br, rect.emission, rect.color);
	}
}






template<typename T, size_t N>
size_t arraySize(T(&)[N]) { return N; }

/******************************************************************
* Check for closest intersection of a ray with the scene;
* Returns true if intersection is found, as well as ray parameter
* of intersection and id of intersected object
*******************************************************************/
bool Intersect_Scene(const Ray &ray, double *t, int *id, Vector *normal) 
{
    const int n = int(sizeof(recs) / sizeof(Rectangle));
    *t  = 1e20;
    *id = -1;
	
    for (int i = 0; i < n; i ++)
    {
        double d = recs[i].intersect(ray);
        if (d > 0.0 && d < *t) 
        {
            *t  = d;
            *id = i;
            *normal = recs[i].normal;
        }
    }
    return *t < 1e20;
}

bool Intersect_Scene_Triangle(const Ray &ray, double *t, int *id, Vector *normal)
{
	const int n = triangles.size();
	bool hit = false;
	*t = FLT_MAX;
	*id = -1;

	for (int i = 0; i < n; i++)
	{
		float d, u, v;
		if (!triangles[i].intersects(ray, d))
			continue;
		else
		{
			hit = true;
			if (d > 0.0 && d < *t)
			{
				*t = d;
				*id = i;
				*normal = triangles[i].getNormal();
			}
		}
	}
	return hit;
}

/******************************************************************
* Determine all form factors for all pairs of patches (of all
* rectangles);
* Evaluation of integrals in form factor equation is done via
* Monte Carlo integration; samples are uniformly distributed and
* equally weighted;
* Computation accelerated by exploiting symmetries of form factor
* estimation;
*******************************************************************/
void Calculate_Form_Factors_Rectangle(const int div_num, const int mc_sample) 
{
	auto a_div_num = div_num;
	auto b_div_num = div_num;

    /* Total number of patches in scene */
    const int n = int(sizeof(recs) / sizeof(Rectangle));
    for (int i = 0; i < n; i ++) 
    {
        recs[i].init_patchs(a_div_num, b_div_num); 
        patch_num += recs[i].a_num * recs[i].b_num;
    }
    
    std::cout << "Number of rectangles: " << n << endl;
    cout << "Number of patches: " << patch_num << endl;
    int form_factor_num = patch_num * patch_num;
    cout << "Number of form factors: " << form_factor_num << endl;
    
    /* 1D-array to hold form factor pairs */
    form_factor = new double[form_factor_num];
    memset(form_factor, 0, sizeof(double) * form_factor_num);

    /* 1D-array with patch areas */
    double *patch_area = new double[patch_num];
    memset(patch_area, 0, sizeof(double) * patch_num);

    /* Precompute patch areas, assuming same size for each rectangle */
    for (int i = 0; i < n; i ++) 
    {
        int patch_i = 0;
        for (int k = 0; k < i; k ++)
            patch_i += recs[k].a_num * recs[k].b_num;

        for (int ia = 0; ia < recs[i].a_num; ia ++) 
        {
            for (int ib = 0; ib < recs[i].b_num; ib ++) 
            {
                patch_area[patch_i + ia* recs[i].b_num + ib] =
					((recs[i].edge_a / recs[i].a_num).
                            Cross((recs[i].edge_b / recs[i].b_num))).Length();  
            }
        }
    }

    /* Offsets for indexing of patches in 1D-array */
    int *offset = new int[n];
    
    for (int i = 0; i < n; i ++) 
    {
        offset[i] = 0;
        for (int k = 0; k < i; k ++)
            offset[i] += recs[k].a_num * recs[k].b_num;
    }

    /* Loop over all rectangles in scene */
    for (int i = 0; i < n; i ++) 
    {
        int patch_i = offset[i];	

        cout << i << " ";
		
        /* Loop over all patches in rectangle i */
        for (int ia = 0; ia < recs[i].a_num; ia ++) 
        {
            cout << "*" << flush;
            for (int ib = 0; ib < recs[i].b_num; ib ++) 
            {
                const Vector normal_i = recs[i].normal;
 
                int patch_j = 0;

                /* Loop over all rectangles in scene for rectangle i */
                for (int j = 0; j < n; j ++) 
                {
                    const Vector normal_j = recs[j].normal;

                    /* Loop over all patches in rectangle j */
                    for (int ja = 0; ja < recs[j].a_num; ja ++) 
                    {
                        for (int jb = 0; jb < recs[j].b_num; jb ++) 
                        {				
                            /* Do not compute form factors for patches on same rectangle;
                               also exploit symmetry to reduce computation;
                               intemediate values; will be divided by patch area below */
                            if (i < j) 
                            {
                                double F = 0;

                                /* Monte Carlo integration of form factor double integral */
                                const int Ni = mc_sample, Nj = mc_sample;
 
                                /* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
                                const double pdf = 
                                    (1.0 / patch_area[offset[i] + ia*recs[i].b_num + ib]) *
                                    (1.0 / patch_area[offset[j] + ja*recs[j].b_num + jb]);

                                /* Determine rays of NixNi uniform samples of patch 
                                   on i to NjxNj uniform samples of patch on j */
                                for (int ias = 0; ias < Ni; ias ++) 
                                {
                                    for (int ibs = 0; ibs < Ni; ibs ++) 
                                    {
                                        for (int jas = 0; jas < Nj; jas ++) 
                                        {
                                            for (int jbs = 0; jbs < Nj; jbs ++) 
                                            {
                                                /* Determine sample points xi, xj on both patches */
                                                const double u0 = (double)(ias + 0.5) / Ni, 
                                                             u1 = (double)(ibs + 0.5) / Ni;
                                                const double u2 = (double)(jas + 0.5) / Nj, 
                                                             u3 = (double)(jbs + 0.5) / Nj;

                                                const Vector xi = recs[i].p0 + 
                                                    recs[i].edge_a * ((double)(ia + u0) / recs[i].a_num) + 
                                                    recs[i].edge_b * ((double)(ib + u1) / recs[i].b_num);
                                                const Vector xj = recs[j].p0 + 
                                                    recs[j].edge_a * ((double)(ja + u2) / recs[j].a_num) +
                                                    recs[j].edge_b * ((double)(jb + u3) / recs[j].b_num);

                                                /* Check for visibility between sample points */
                                                const Vector ij = (xj - xi).Normalized();

                                                double t; 
                                                int id;
                                                Vector normal; 
                                                if (Intersect_Scene(Ray(xi, ij), &t, &id, &normal) && 
                                                    id != j) 
                                                {
                                                    continue; /* If intersection with other rectangle */
                                                }

                                                /* Cosines of angles beteen normals and ray inbetween */
                                                const double d0 = normal_i.Dot(ij);
                                                const double d1 = normal_j.Dot(-1.0 * ij);

                                                /* Continue if patches facing each other */
                                                if (d0 > 0.0 && d1 > 0.0) 
                                                {
                                                    /* Sample form factor */
                                                    const double K = d0 * d1 / 
                                                            (M_PI * (xj - xi).LengthSquared());

                                                    /* Add weighted sample to estimate */
                                                    F += K / pdf;
                                                }
                                            }
                                        }
                                    }
                                } 

                                /* Divide by number of samples */
                                F /= (Ni) * (Ni) * (Nj) * (Nj); 
 
                                form_factor[patch_i * patch_num + patch_j] = F;
                            }
                            patch_j ++;
                        }
                    }
                }
                patch_i ++;
            }
        }

        cout << endl;
    }

    /* Copy upper to lower triangular values */
    for (int i = 0; i < patch_num-1; i ++) 
    {
        for (int j = i+1; j < patch_num; j ++) 
        {
             form_factor[j * patch_num + i] = form_factor[i * patch_num + j];
        }
    }

    /* Divide by area to get final form factors */
    for (int i = 0; i < patch_num; i ++) 
    {
        for (int j = 0; j < patch_num; j ++) 
        {
             form_factor[i * patch_num + j] /= patch_area[i];

             /* Clamp to [0,1] */
             if(form_factor[i * patch_num + j] > 1.0) 
                 form_factor[i * patch_num + j] = 1.0;
        }
    }
}

void Calculate_Form_Factors_Triangle(const int div_num, const int mc_sample)
{
	const auto N = triangles.size();

	// Total number of patches in scene
	// each triangle is split into four subtriangles per division
	patch_num = (int)N * div_num * 4; 

	std::cout << "Number of triangles: " << N << endl;
	cout << "Number of patches: " << patch_num << endl;

	auto form_factor_num = patch_num * patch_num;
	cout << "Number of form factors: " << form_factor_num << endl;

	/* 1D-array to hold form factor pairs */
	form_factor = new double[form_factor_num];
	memset(form_factor, 0, sizeof(double) * form_factor_num);

	/* 1D-array with patch areas */
	auto patch_area = new double[patch_num];
	memset(patch_area, 0, sizeof(double) * patch_num);

	/* Precompute patch areas, assuming same size for each rectangle */
	auto subTriangleCount = div_num * 4;
	for (auto i = 0u; i < N; i++)
	{
		auto patch_i = i * subTriangleCount;

		for (auto subTriangle = 0; subTriangle < subTriangleCount; subTriangle++)
		{
			patch_area[patch_i + subTriangle] = triangles[i].area() / subTriangleCount;
		}
	}

	/* Offsets for indexing of patches in 1D-array */
	auto offset = new int[N];

	/* Loop over all rectangles in scene */
	for (auto i = 0u; i < N; i++)
	{
		const auto patch_i = i * subTriangleCount;
		const Vector normal_i = triangles[i].normal();

		cout << i << " ";

		auto sub_triangles_i = new Triangle[subTriangleCount];
		triangles[i].split_n(div_num, sub_triangles_i);

		/* Loop over all patches in rectangle i */
		for (auto sub_triangle_i = 0; sub_triangle_i < subTriangleCount; sub_triangle_i++)
		{
			/* Loop over all rectangles in scene for rectangle i */
			for (auto j = 0u; j < N; j++)
			{
				const auto patch_j = j * subTriangleCount;
				const Vector normal_j = triangles[j].normal();

				auto sub_triangles_j = new Triangle[subTriangleCount];
				triangles[j].split_n(div_num, sub_triangles_j);
				
				/* Loop over all patches in rectangle j */
				for (auto sub_triangle_j = 0; sub_triangle_j < subTriangleCount; sub_triangle_j++)
				{
					/* Do not compute form factors for patches on same rectangle;
					also exploit symmetry to reduce computation;
					intemediate values; will be divided by patch area below */
					if (i < j)
					{
						double F = 0;

						/* Monte Carlo integration of form factor double integral */
						const int Ni = mc_sample, Nj = mc_sample;

						/* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
						const double pdf =
							(1.0 / patch_area[patch_i + sub_triangle_i]) *
							(1.0 / patch_area[patch_j + sub_triangle_j]);

						///////////////////////////////////////////////////////////////
						// simple form factor calculcation with only one sample
						const auto xi = sub_triangles_i[sub_triangle_i].point_inside();
						const auto xj = sub_triangles_j[sub_triangle_j].point_inside();

						// Check for visibility between sample points
						const auto ij = (xj - xi).Normalized();

						double t;
						int id;
						Vector normal;
						if (Intersect_Scene_Triangle(Ray(xi, ij), &t, &id, &normal) &&
							id != j)
						{
							continue; //If intersection with other rectangle
						}

						// Cosines of angles beteen normals and ray inbetween
						const double d0 = normal_i.Dot(ij);
						const double d1 = normal_j.Dot(-1.0 * ij);

						// Continue if patches facing each other
						if (d0 > 0.0 && d1 > 0.0)
						{
							// Sample form facto´r
							const double K = d0 * d1 /
								(M_PI * (xj - xi).LengthSquared());

							//Add weighted sample to estimate
							F += K / pdf;
						}

						///////////////////////////////////////////////////////////////

						/* Determine rays of NixNi uniform samples of patch
						on i to NjxNj uniform samples of patch on j
						/*for (int ias = 0; ias < Ni; ias++)
						{
							for (int ibs = 0; ibs < Ni; ibs++)
							{
								for (int jas = 0; jas < Nj; jas++)
								{
									for (int jbs = 0; jbs < Nj; jbs++)
									{
										// Determine sample points xi, xj on both patche
										const auto u0 = (ias + 0.5) / Ni;
										const auto u1 = (ibs + 0.5) / Ni;
										const auto u2 = (jas + 0.5) / Nj;
										const auto u3 = (jbs + 0.5) / Nj;

										const Vector xi = recs[i].p0 +
											recs[i].edge_a * ((double)(ia + u0) / recs[i].a_num) +
											recs[i].edge_b * ((double)(ib + u1) / recs[i].b_num);

										const Vector xj = recs[j].p0 +
											recs[j].edge_a * ((double)(ja + u2) / recs[j].a_num) +
											recs[j].edge_b * ((double)(jb + u3) / recs[j].b_num);

										// Check for visibility between sample points
										const auto ij = (xj - xi).Normalized();

										double t;
										int id;
										Vector normal;
										if (Intersect_Scene_Triangle(triangles, Ray(xi, ij), &t, &id, &normal) &&
											id != j)
										{
											continue; //If intersection with other rectangle
										}

										// Cosines of angles beteen normals and ray inbetween
										const double d0 = normal_i.Dot(ij);
										const double d1 = normal_j.Dot(-1.0 * ij);

										// Continue if patches facing each other
										if (d0 > 0.0 && d1 > 0.0)
										{
											// Sample form facto´r
											const double K = d0 * d1 /
												(M_PI * (xj - xi).LengthSquared());

											//Add weighted sample to estimate
											F += K / pdf;
										}
									}
								}
							}
						}

						// Divide by number of samples
						F /= (Ni)* (Ni)* (Nj)* (Nj);
						
						*/

						form_factor[patch_i * patch_num + patch_j] = F;
					}
				}

				delete[] sub_triangles_j;
			}
		}

		delete[] sub_triangles_i;

		cout << endl;
	}

	/* Copy upper to lower triangular values */
	for (auto i = 0; i < patch_num - 1; i++)
	{
		for (auto j = i + 1; j < patch_num; j++)
		{
			form_factor[j * patch_num + i] = form_factor[i * patch_num + j];
		}
	}

	/* Divide by area to get final form factors */
	for (auto i = 0; i < patch_num; i++)
	{
		for (auto j = 0; j < patch_num; j++)
		{
			form_factor[i * patch_num + j] /= patch_area[i];

			/* Clamp to [0,1] */
			if (form_factor[i * patch_num + j] > 1.0)
				form_factor[i * patch_num + j] = 1.0;
		}
	}

	delete[] offset;
}

/******************************************************************
* Iterative computation of radiosity via Gathering; i.e. solution
* using Gauss-Seidel iteration - reuse already computed values;
* run-time O(n^2) 
*******************************************************************/

void Calculate_Radiosity_Rectangle(const int iteration)
{
    const int n = int(sizeof(recs) / sizeof(Rectangle));
    int patch_i = 0;
	
    for (int i = 0; i < n; i ++) 
    {
        for (int ia = 0; ia < recs[i].a_num; ia ++) 
        {
            for (int ib = 0; ib < recs[i].b_num; ib ++) 
            {
                Color B;

                int patch_j = 0;
                for (int j = 0; j < n; j ++) 
                {
                    for (int ja = 0; ja < recs[j].a_num; ja ++) 
                    {
                        for (int jb = 0; jb < recs[j].b_num; jb ++) 
                        {
                            const double Fij = form_factor[patch_i * patch_num + patch_j];
					 
                            /* Add form factor multiplied with radiosity of previous step */		
                            if (Fij > 0.0)
                                B = B + Fij * recs[j].patch[ja * recs[j].b_num + jb];

                            patch_j ++;
                        }
                    }
                }
                /* Multiply sum with color of patch and add emission */
                B = recs[i].color.MultComponents(B) + recs[i].emission;
  
                /* Store overall patch radiosity of current iteration */
                recs[i].patch[ia * recs[i].b_num + ib] = B;
                patch_i ++;
            }
        }
    }
}

void Calculate_Radiosity_Triangle(const int iteration) { }

/******************************************************************
* Helper functions for smooth bicubic (Catmull-Rom) interpolation 
* using 4x4 color patches;
* First interpolate in y, followed by interpolation of results in x
*******************************************************************/

Color cubicInterpolate (Color p[4], double x) 
{
    return p[1] + 0.5 * x * (p[2] - p[0] + x * (2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + 
                                                 x * (3.0*(p[1] - p[2]) + p[3] - p[0])));
}

Color bicubicInterpolate (Color p[4][4], double x, double y) 
{
	Color arr[4];

	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	
    return cubicInterpolate(arr, x);
}


/******************************************************************
* Compute radiance from radiosity by shooting rays into the scene;
* Radiance directly proportional to radiosity for assumed diffuse
* emitters/surfaces (multiply by PI);
* At intersections either constant patch color is returned or a
* smoothly interpolated color of 4x4 neighboring patches
*******************************************************************/

Color Radiance_Rectangle(const Ray &ray, const int depth, bool interpolation = true)
{
    double t; 
    int id;  
    Vector normal; 

    /* Find intersected rectangle */
    if (!Intersect_Scene(ray, &t, &id, &normal))
    {
        return BackgroundColor;    
    }

    /* Determine intersection point on rectangle */
    const Rectangle &obj = recs[id];
    const Vector hitpoint = ray.org + t * ray.dir; 

    /* Determine intersected patch */
    const Vector v = hitpoint - obj.p0;
    const double a_len = v.Dot(obj.edge_a.Normalized());
    const double b_len = v.Dot(obj.edge_b.Normalized());
            
    double da = obj.a_num * a_len / obj.a_len;
    double db = obj.b_num * b_len / obj.b_len;

    int ia = int(da); if (ia >= obj.a_num) ia --;
    int ib = int(db); if (ib >= obj.b_num) ib --;
            
    /* Bicubic interpolation for smooth image */
    if (interpolation)  
    {
        Color c[4][4];

        int ia = int(da - 0.5);
        int ib = int(db - 0.5);
                
        for (int i = 0; i < 4; i ++) 
        {
            for (int j = 0; j < 4; j ++) 
            {
                c[i][j] = obj.sample_patch(ia + i - 1, ib + j - 1);
            }
        }
               
        int ia0 = int(da - 0.5);
        int ib0 = int(db - 0.5); 
        double dx = da - ia0 - 0.5;
        double dy = db - ib0 - 0.5;

        if (dx < 0.0)  dx = 0.0;
        if (dx >= 1.0) dx = 1.0;
        if (dy < 0.0)  dy = 0.0;
        if (dy >= 1.0) dy = 1.0;
 
        return bicubicInterpolate(c, dx, dy) * Over_M_PI;
    }
    else
    {         
        return obj.patch[ia * obj.b_num + ib] * Over_M_PI;
    }
}

Color Radiance_triangles(const Ray &ray, const int depth, bool interpolation = true)
{
	double t;
	int id;
	Vector normal;

	/* Find intersected rectangle */
	if (!Intersect_Scene_Triangle(ray, &t, &id, &normal))
	{
		return BackgroundColor;
	}

	/* Determine intersection point on rectangle */
	const Triangle &obj = triangles[id];
	const Vector hitpoint = ray.org + t * ray.dir;


	/* Determine intersected patch */
	const Vector v = hitpoint - obj.getP1();
	auto edgeA = obj.getP2() - obj.getP1();
	auto edgeB = obj.getP3() - obj.getP1();

	const double a_len = v.Dot(edgeA.Normalized());
	const double b_len = v.Dot(edgeB.Normalized());

	return obj.getColor();

// 	double da = obj.a_num * a_len / obj.a_len;
// 	double db = obj.b_num * b_len / obj.b_len;
// 
// 	int ia = int(da); if (ia >= obj.a_num) ia--;
// 	int ib = int(db); if (ib >= obj.b_num) ib--;
// 
// 	/* Bicubic interpolation for smooth image */
// 	if (interpolation)
// 	{
// 		Color c[4][4];
// 
// 		int ia = int(da - 0.5);
// 		int ib = int(db - 0.5);
// 
// 		for (int i = 0; i < 4; i++)
// 		{
// 			for (int j = 0; j < 4; j++)
// 			{
// 				c[i][j] = obj.sample_patch(ia + i - 1, ib + j - 1);
// 			}
// 		}
// 
// 		int ia0 = int(da - 0.5);
// 		int ib0 = int(db - 0.5);
// 		double dx = da - ia0 - 0.5;
// 		double dy = db - ib0 - 0.5;
// 
// 		if (dx < 0.0)  dx = 0.0;
// 		if (dx >= 1.0) dx = 1.0;
// 		if (dy < 0.0)  dy = 0.0;
// 		if (dy >= 1.0) dy = 1.0;
// 
// 		return bicubicInterpolate(c, dx, dy) * Over_M_PI;
// 	}
// 	else
// 	{
// 		return obj.patch[ia * obj.b_num + ib] * Over_M_PI;
// 	}
}

#define USE_TRIANGLES

#ifdef USE_TRIANGLES
#define Calculate_Form_Factors Calculate_Form_Factors_Triangle
#define Calculate_Radiosity Calculate_Radiosity_Triangle
#define Radiance Radiance_triangles
#else
#define Calculate_Form_Factors Calculate_Form_Factors_Rectangle
#define Calculate_Radiosity Calculate_Radiosity_Rectangle
#define Calculate_Radiosity Calculate_Radiosity_Rectangle
#define Radiance Radiance_Rectangle
#endif

/******************************************************************
* Main routine: Computation of radiosity image
* Key parameters
* - Image dimensions: width, height 
* - Number of samples for antialiasing (non-uniform filter): samples 
* - Number of patches along edges a,b: patches_a, patches_b
* - Number of uniform samples per patch edge: MC_samples
* - Number of iterations for iterative solver: iterations
* Rendered result saved as PPM image file
*******************************************************************/

int main(int argc, char **argv) 
{
	convertRectangles(recs, triangles);

    int width = 640;
    int height = 480;
    int samples = 2;

    /* Set camera origin and viewing direction (negative z direction) */
    Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());

    /* Image edge vectors for pixel sampling */
    Vector cx = Vector(width * 0.5135 / height);
    Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

    /* Two final renderings; one with constant, one with interpolated patch colors */
    Image img(width, height);
    Image img_interpolated(width, height);

    cout << "Calculating form factors" << endl;
    int divisions = 4;
    int MC_samples = 3;

    Calculate_Form_Factors(divisions, MC_samples);

    /* Iterative solution of radiosity linear system */
    cout << "Calculating radiosity" << endl;
    int iterations = 40; 
    for (int i = 0; i < iterations; i ++) 
    {
        cout << i << " ";
        Calculate_Radiosity(i);
    }
    cout << endl;
 
    /* Loop over image rows */
    for (int y = 0; y < height; y ++) 
    {
        cout << "\rRendering (" << samples * 4 << " spp) " << (100.0 * y / (height - 1)) << "%     ";
        srand(y * y * y);

        /* Loop over row pixels */
        for (int x = 0; x < width; x ++) 
        {
            img.setColor(x, y, Color());
            img_interpolated.setColor(x, y, Color());

            /* 2x2 subsampling per pixel */
            for (int sy = 0; sy < 2; sy ++) 
            {
                for (int sx = 0; sx < 2; sx ++) 
                {
                    Color accumulated_radiance = Color();
                    Color accumulated_radiance2 = Color();

                    /* Computes radiance at subpixel using multiple samples */
                    for (int s = 0; s < samples; s ++) 
                    {
                        const double r1 = 2.0 * drand48();
                        const double r2 = 2.0 * drand48();

                        /* Transform uniform into non-uniform filter samples */
                        double dx;
                        if (r1 < 1.0)
                            dx = sqrt(r1) - 1.0;
                        else
                            dx = 1.0 - sqrt(2.0 - r1);

                        double dy;
                        if (r2 < 1.0)
                            dy = sqrt(r2) - 1.0;
                        else
                            dy = 1.0 - sqrt(2.0 - r2);

                        /* Ray direction into scene from camera through sample */
                        Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / width - 0.5) +
                                     cy * ((y + (sy + 0.5 + dy) / 2.0) / height - 0.5) + 
                                     camera.dir;
                        
                        /* Extend camera ray to start inside box */
                        Vector start = camera.org + dir * 130.0;

                        /* Determine constant radiance */
                        accumulated_radiance = accumulated_radiance + 
                            Radiance(Ray(start, dir.Normalized()), 0, false) / samples;

                        /* Determine interpolated radiance */
                        accumulated_radiance2 = accumulated_radiance2 + 
                            Radiance(Ray(start, dir.Normalized()), 0, true) / samples;
                    }

                    img.addColor(x, y, accumulated_radiance);
                    img_interpolated.addColor(x, y, accumulated_radiance2); 
                 }
            }
        }
    }

    cout << endl;
	
    img.Save(string("image_patches.ppm"));
    img_interpolated.Save(string("image_smooth.ppm"));

#ifdef _WIN32
	system("image_smooth.ppm");
#endif
}
