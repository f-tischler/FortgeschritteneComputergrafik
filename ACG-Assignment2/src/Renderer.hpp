#ifndef Renderer_H_included
#define Renderer_H_included

#include <memory>
#include <vector>

#include "Consts.hpp"
#include "Geometry.hpp"
#include "Sphere.hpp"
#include "TriangleMesh.hpp"
#include "Image.hpp"
#include "Scene.hpp"


class Renderer
{
protected:
    unsigned int mWidth;
    unsigned int mHeight;
    unsigned int mSamples;
    Color mClearColor;
    SPScene mScene;
    
    //Safe scope
    IntersectionInfo mIntersectionInfo;
    
public:
    Renderer() : mClearColor(0, 1, 1)
    {
    }
    
    void buildScene(SPScene scene)
    {
        mScene = scene;
    }
    
    /******************************************************************
     * Check for closest intersection of a ray with the scene;
     * returns true if intersection is found, as well as ray parameter
     * of intersection and id of intersected object
     *******************************************************************/
    bool Intersect(const Ray &ray, double &t, size_t& id, Vector& normal, Vector& hitpoint, bool culling = true)
    {
        Vector normalIntersect;
        t = 1e20;
        if (!mScene->Intersect(ray, mIntersectionInfo, culling))
            return false;
        
        normal = mIntersectionInfo.normal;
        hitpoint = mIntersectionInfo.hitpoint;
        id = mIntersectionInfo.geometryId;
        t = mIntersectionInfo.distance;
        
        return t < 1e20;
    }
    
    Image render(unsigned int width, unsigned int height, unsigned int samples, Color clearColor = Color(0, 1, 1))
    {
        mClearColor = clearColor;
        
        mWidth = width;
        mHeight = height;
        mSamples = samples;
        
        /* Set camera origin and viewing direction (negative z direction) */
        Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());
        
        /* Image edge vectors for pixel sampling */
        Vector cx = Vector(mWidth * 0.5135 / mHeight);
        Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;
        
        /* Final rendering */
        Image img(mWidth, mHeight);
        
        /* Loop over image rows */
        for (auto y = 0u; y < mHeight; y++)
        {
            std::cout << "\rRendering (" << mSamples * 4 << " spp) " << (100.0 * y / (mHeight - 1)) << "%     ";
            srand(y * y * y);
            
            /* Loop over row pixels */
            for (auto x = 0u; x < mWidth; x++)
            {
                img.setColor(x, y, Color());
                
                /* 2x2 subsampling per pixel */
                for (int sy = 0; sy < 2; sy++)
                {
                    for (int sx = 0; sx < 2; sx++)
                    {
                        Color accumulated_radiance = Color();
                        
                        /* Compute radiance at subpixel using multiple samples */
                        for (auto s = 0u; s < mSamples; s++)
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
                            Vector dir = cx * ((x + (sx + 0.5 + dx) / 2.0) / mWidth - 0.5) +
                            cy * ((y + (sy + 0.5 + dy) / 2.0) / mHeight - 0.5) +
                            camera.dir;
                            
                            /* Extend camera ray to start inside box */
                            Vector start = camera.org + dir * 130.0;
                            
                            dir = dir.Normalized();
                            
                            accumulated_radiance = accumulated_radiance + Radiance(Ray(start, dir), 0, 1) / mSamples;
                        }
                        
                        accumulated_radiance = accumulated_radiance.clamp() * 0.25;
                        
                        img.addColor(x, y, accumulated_radiance);
                    }
                }
            }
        }
        std::cout << std::endl;
        
        return img;
    }
    
private:
    /******************************************************************
     * Recursive path tracing for computing radiance via Monte-Carlo
     * integration; only considers perfectly diffuse, specular or
     * transparent materials;
     * after 5 bounces Russian Roulette is used to possibly terminate rays;
     * emitted light from light source only included on first direct hit
     * (possibly via specular reflection, refraction), controlled by
     * parameter E = 0/1;
     * on diffuse surfaces light sources are explicitely sampled;
     * for transparent objects, Schlick´s approximation is employed;
     * for first 3 bounces obtain reflected and refracted component,
     * afterwards one of the two is chosen randomly
     *******************************************************************/
    Color Radiance(const Ray &ray, int depth, int E)
    {
        depth++;
        
        if (!mScene->Intersect(ray, mIntersectionInfo))   /* No intersection with scene */
            return mClearColor;
        
        //const auto &obj = *(mIntersectionInfo.geometry.get());
        const auto &obj = *mIntersectionInfo.geometry;
        
        Vector nl = mIntersectionInfo.normal;

			
        /* Obtain flipped normal, if object hit from inside */
        if (mIntersectionInfo.normal.Dot(ray.dir) >= 0)
            nl = nl*-1.0;
        
        Color col = obj.color;
        
        /* Maximum RGB reflectivity for Russian Roulette */
        double p = col.Max();
        
        if (depth > 5 || !p)   /* After 5 bounces or if max reflectivity is zero */
        {
            if (drand48() < p*0.9)            /* Russian Roulette */
                col = col * (1 / p);        /* Scale estimator to remain unbiased */
            else
                return obj.emission * E;  /* No further bounces, only return potential emission */
        }
        
        if (obj.refl == DIFF)
        {
            /* Compute random reflection vector on hemisphere */
            double r1 = 2.0 * M_PI * drand48();
            double r2 = drand48();
            double r2s = sqrt(r2);
            
            /* Set up local orthogonal coordinate system u,v,w on surface */
            Vector w = nl;
            Vector u;
            
            if (fabs(w.x) > .1)
                u = Vector(0.0, 1.0, 0.0);
            else
                u = (Vector(1.0, 0.0, 0.0).Cross(w)).Normalized();
            
            Vector v = w.Cross(u);
            
            /* Random reflection vector d */
            Vector d = (u * cos(r1) * r2s +
                        v * sin(r1) * r2s +
                        w * sqrt(1 - r2)).Normalized();
            
            /* Explicit computation of direct lighting */
            Vector e;
            for (auto i = 0u; i < mScene->geometries().size(); i++)
            {
                auto geom = mScene->geometries()[i];
                if (geom->emission.x <= 0 && geom->emission.y <= 0 && geom->emission.z <= 0)
                    continue;
                
                if (geom->getType() == eGeometryType_Sphere)
                {
                    const Sphere &sphere = *(Sphere*)(geom.get());
                    /* Randomly sample spherical light source from surface intersection */
                    
                    /* Set up local orthogonal coordinate system su,sv,sw towards light source */
                    Vector sw = sphere.position - mIntersectionInfo.hitpoint;
                    Vector su;
                    
                    if (fabs(sw.x) > 0.1)
                        su = Vector(0.0, 1.0, 0.0);
                    else
                        su = Vector(1.0, 0.0, 0.0);
                    
                    su = (su.Cross(w)).Normalized();
                    Vector sv = sw.Cross(su);
                    
                    /* Create random sample direction l towards spherical light source */
                    double cos_a_max = sqrt(1.0 - sphere.radius * sphere.radius /
                                            (mIntersectionInfo.hitpoint - sphere.position).Dot(mIntersectionInfo.hitpoint - sphere.position));
                    double eps1 = drand48();
                    double eps2 = drand48();
                    double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
                    double sin_a = sqrt(1.0 - cos_a * cos_a);
                    double phi = 2.0*M_PI * eps2;
                    Vector l = su * cos(phi) * sin_a +
                    sv * sin(phi) * sin_a +
                    sw * cos_a;
                    l = l.Normalized();
                    
                    /* Shoot shadow ray, check if intersection is with light source */
                    IntersectionInfo info;
                    if (mScene->Intersect(Ray(mIntersectionInfo.hitpoint, l), info) && info.geometryId == i)
                    {
                        double omega = 2 * M_PI * (1 - cos_a_max);
                        
                        /* Add diffusely reflected light from light source; note constant BRDF 1/PI */
                        e = e + col.MultComponents(sphere.emission * l.Dot(nl) * omega) * M_1_PI;
                    }
                }
                else
                {
                    const auto &mesh = *(TriangleMesh*)(geom.get());
                    /* Randomly sample spherical light source from surface intersection */
                    auto radius = mesh.getBB().diameter()/2;
                    auto position = mesh.getBB().position();
                    
                    /* Set up local orthogonal coordinate system su,sv,sw towards light source */
                    Vector sw = position - mIntersectionInfo.hitpoint;
                    Vector su;
                    
                    if (fabs(sw.x) > 0.1)
                        su = Vector(0.0, 1.0, 0.0);
                    else
                        su = Vector(1.0, 0.0, 0.0);
                    
                    su = (su.Cross(w)).Normalized();
                    Vector sv = sw.Cross(su);
                    
                    /* Create random sample direction l towards spherical light source */
                    double cos_a_max = sqrt(1.0 - radius * radius /
                                            (mIntersectionInfo.hitpoint - position).Dot(mIntersectionInfo.hitpoint - position));
                    double eps1 = drand48();
                    double eps2 = drand48();
                    double cos_a = 1.0 - eps1 + eps1 * cos_a_max;
                    double sin_a = sqrt(1.0 - cos_a * cos_a);
                    double phi = 2.0*M_PI * eps2;
                    Vector l = su * cos(phi) * sin_a +
                    sv * sin(phi) * sin_a +
                    sw * cos_a;
                    l = l.Normalized();
                    
                    /* Shoot shadow ray, check if intersection is with light source */
                    IntersectionInfo info;
                    if (mScene->Intersect(Ray(mIntersectionInfo.hitpoint, l), info) && info.geometryId == i)
                    {
                        double omega = 2 * M_PI * (1 - cos_a_max);
                        
                        /* Add diffusely reflected light from light source; note constant BRDF 1/PI */
                        e = e + col.MultComponents(mesh.emission * l.Dot(nl) * omega) * M_1_PI;
                    }
                }
            }
            
            /* Return potential light emission, direct lighting, and indirect lighting (via
             recursive call for Monte-Carlo integration */
            return obj.emission * E + e + col.MultComponents(Radiance(Ray(mIntersectionInfo.hitpoint, d), depth, 0));
        }
        else if (obj.refl == SPEC)
        {
            /* Return light emission mirror reflection (via recursive call using perfect 
             reflection vector) */
            return col.MultComponents(Radiance(Ray(mIntersectionInfo.hitpoint, ray.dir - mIntersectionInfo.normal * 2 * mIntersectionInfo.normal.Dot(ray.dir)),depth, 1));
        }
        else if (obj.refl == GLOS)
        {
            Vector L = ray.dir;
            Vector N = mIntersectionInfo.normal;
            Vector L_prime = L - N * 2 * (N.Dot(L));
         
            varyVector(L_prime, obj.getConstant());
         
            return obj.emission + col.MultComponents(Radiance(Ray(mIntersectionInfo.hitpoint, L_prime), depth, 1));
        }
        
		/* Otherwise object transparent, i.e. assumed dielectric glass material */
		Vector L = ray.dir;
		Vector normal = mIntersectionInfo.normal;
		auto hitpoint = mIntersectionInfo.hitpoint;
		Vector N = normal;
		Vector L_prime = L - N * 2 * (N.Dot(L));
		if (obj.refl == TRAN) varyVector(L_prime, 0.8);
		Ray reflRay(hitpoint, L_prime);  /* Prefect reflection */
		bool into = normal.Dot(nl) > 0;       /* Bool for checking if ray from outside going in */
		double nc = 1;                        /* Index of refraction of air (approximately) */
		double nt = 1.5;                      /* Index of refraction of glass (approximately) */
		double nnt;

		if (into)      /* Set ratio depending on hit from inside or outside */
			nnt = nc / nt;
		else
			nnt = nt / nc;

		double ddn = ray.dir.Dot(nl);
		double cos2t = 1 - nnt * nnt * (1 - ddn*ddn);

		/* Check for total internal reflection, if so only reflect */
		if (cos2t < 0)
			return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1));

		/* Otherwise reflection and/or refraction occurs */
		Vector tdir;

		/* Determine transmitted ray direction for refraction */
		if (into)
			tdir = (ray.dir * nnt - normal * (ddn * nnt + sqrt(cos2t))).Normalized();
		else
			tdir = (ray.dir * nnt + normal * (ddn * nnt + sqrt(cos2t))).Normalized();


		if (obj.refl == TRAN) varyVector(tdir, 0.8);

		/* Determine R0 for Schlick´s approximation */
		double a = nt - nc;
		double b = nt + nc;
		double R0 = a*a / (b*b);

		/* Cosine of correct angle depending on outside/inside */
		double c;
		if (into)
			c = 1 + ddn;
		else
			c = 1 - tdir.Dot(normal);

		/* Compute Schlick´s approximation of Fresnel equation */
		double Re = R0 + (1 - R0) *c*c*c*c*c;   /* Reflectance */
		double Tr = 1 - Re;                     /* Transmittance */

		/* Probability for selecting reflectance or transmittance */
		double P = .25 + .5 * Re;
		double RP = Re / P;         /* Scaling factors for unbiased estimator */
		double TP = Tr / (1 - P);

		if (depth < 3)   /* Initially both reflection and trasmission */
			return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * Re +
			Radiance(Ray(hitpoint, tdir), depth, 1) * Tr);
		else             /* Russian Roulette */
		if (drand48() < P)
			return obj.emission + col.MultComponents(Radiance(reflRay, depth, 1) * RP);
		else
			return obj.emission + col.MultComponents(Radiance(Ray(hitpoint, tdir), depth, 1) * TP);
	}
    
    // Cumulative distribution function
    // c can be used as glossiness/translucency factor. c in interval (0, infinity)
    // small c (e.g. 0.01) means mirror like, large c (e.g. 0.9) means really glossy
    // x should be a random number in interval [0, 1]
    double CDF(double x, double c) {
        c = c <= 0 ? 1.0e-6 : c;
        double beta = atan(-1.0/c);
        double alpha = atan(1.0/c) - beta;
        return c * tan(alpha * x + beta);
    }
    

	// myFunction1
	// input:  x should be a random number between 0 and 1
	//         glossiness should be a value between 0 and 1
	//         small glossiness means really blurry reflection (evenly distributed randomness)
	//         large glossiness means really sharp reflection (unevenly distributed randomness)
	// output: a value between -1 and 1
	double myFunction1(double x, double glossiness) {
		if (glossiness < 0) glossiness = 0;
		if (glossiness >= 1.0) glossiness = 1.0 - 1.0e-6;
		double c = -1.0 / log(glossiness);
		double beta = atan(-c);
		double alpha = atan(c) - beta;
		return (1.0 / c) * tan(alpha * x + beta);
	}

	// myFunction2
	// input:  x should be a random number between 0 and 1
	//         glossiness should be a value between 0 and 1
	//         small glossiness means really blurry reflection (evenly distributed randomness)
	//         large glossiness means really sharp reflection (unevenly distributed randomness)
	// output: a value between 0 and 1
	double myFunction2(double x, double glossiness) {
		if (glossiness < 0) glossiness = 0;
		if (glossiness >= 1.0) glossiness = 1.0 - 1.0e-6;
		double c = -1.0 / log(glossiness);
		return exp(c * x - c);
	}

	// varyVector varies vector
	void varyVector(Vector &vector, double glossiness) {
		/* Compute random reflection vector on half hemisphere */
		double oldtheta = acos(vector.z);
		double oldphi = atan(vector.y / vector.x);
		double deltatheta = M_PI * (2 * drand48() - 1);
		double deltaphi = (M_PI / 4.0) * (2 * drand48() - 1);

		double theta = oldtheta + deltatheta;
		double phi = oldphi + deltaphi;

		/* Random reflection vector d */
		Vector d = Vector(sin(theta) * cos(phi),
			sin(theta) * sin(phi),
			cos(theta));

		vector += myFunction2(drand48(), glossiness) * d; // myFunction2 works better than myFunction1, since the vector will still point in the same general direction
		vector = vector.Normalized();
	}
};

#endif