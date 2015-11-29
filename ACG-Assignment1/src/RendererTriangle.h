#ifndef RendererTriangles_H_included
#define RendererTriangles_H_included

#include <vector>

#include "Renderer.hpp"
#include "Rectangle.hpp"
#include "Triangle.hpp"

#define drand48() (((double)rand())/((double)RAND_MAX))

class RendererTriangles : public Renderer
{
	const double Over_M_PI = 1.0 / M_PI;

	std::vector<double> mFormFactor;
	unsigned int mPatchCount;
	std::vector<Triangle> mTriangles;

public:
	RendererTriangles(unsigned int mWidth, unsigned int mHeight, unsigned int mSamples) : Renderer(mWidth, mHeight, mSamples)
	{
		mPatchCount = 0;
	}

	virtual void buildScene(const std::vector<Rectangle>& rects)
	{
		if (!mTriangles.empty())
			mTriangles.clear();

		mTriangles.reserve(rects.size() * 2);
		for (auto rect : rects)
		{
			auto tl = rect.p0;
			auto bl = rect.p0 + rect.edge_b;
			auto tr = rect.p0 + rect.edge_a;
			auto br = bl + rect.edge_a;

			auto t1 = Triangle(tl, tr, bl, rect.emission, rect.color);
			auto t2 = Triangle(bl, tr, br, rect.emission, rect.color);
			
			//std::cout << "t1: " << t1.normal().toString() << std::endl;	
			//std::cout << "t2: " << t2.normal().toString() << std::endl;

			const auto epsilon = 0.00001f;
			const auto epsilonVec = Vector(epsilon, epsilon, epsilon);

			assert(Vector::AreEqual(rect.normal,
									t1.normal(),
									epsilonVec) &&
					"triangle normals must match rectangle normal");

			assert(Vector::AreEqual(t1.normal(), 
									t2.normal(), 
									epsilonVec) &&
				   "normals for the rectangle must match");

			mTriangles.emplace_back(std::move(t1));
			mTriangles.emplace_back(std::move(t2));
		}
	}

	virtual void Render(Image& img, Image& imgInterpolated)
	{
		/* Set camera origin and viewing direction (negative z direction) */
		Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());

		/* Image edge vectors for pixel sampling */
		auto cx = Vector(mWidth * 0.5135 / mHeight);
		auto cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

		std::cout << "Calculating form factors" << std::endl;
		const auto divisions = 1; // subtriangleCount = 4^divisions
		const int MC_mSamples = 3;

		Calculate_Form_Factors(divisions, MC_mSamples);

		/* Iterative solution of radiosity linear system */
		std::cout << "Calculating radiosity" << std::endl;
		const auto iterations = 10;
		for (auto i = 0; i < iterations; i++)
		{
			std::cout << i << " ";
			Calculate_Radiosity(i);
		}
		std::cout << std::endl;

		/* Loop over image rows */
		for (auto y = 0u; y < mHeight; y++)
		{
			std::cout << "\rRendering (" << mSamples * 4 << " spp) " << (100.0 * y / (mHeight - 1)) << "%     ";
			srand(y * y * y);

			/* Loop over row pixels */
			for (auto x = 0u; x < mWidth; x++)
			{
				img.setColor(x, y, Color());
				imgInterpolated.setColor(x, y, Color());

				/* 2x2 subsampling per pixel */
				for (auto sy = 0u; sy < 2; sy++)
				{
					for (int sx = 0; sx < 2; sx++)
					{
						auto accumulated_radiance = Color();
						auto accumulated_radiance2 = Color();

						/* Computes radiance at subpixel using multiple mSamples */
						for (auto s = 0u; s < mSamples; s++)
						{
							const double r1 = 2.0 * drand48();
							const double r2 = 2.0 * drand48();

							/* Transform uniform into non-uniform filter mSamples */
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
							auto ray = Ray(start, dir.Normalized());

							/* Find intersected rectangle */
							double distance;
							int id;
							Vector normal;
							if (!Intersect_Scene(ray, &distance, &id, &normal))
							{
								accumulated_radiance += mClearColor;
								accumulated_radiance2 += mClearColor;
							}
							else
							{
								/* Determine constant radiance */
								accumulated_radiance += Radiance(id, distance, ray, 0, false) / mSamples;

								/* Determine interpolated radiance */
								accumulated_radiance2 += Radiance(id, distance, ray, 0, true) / mSamples;
							}
						}

						img.addColor(x, y, accumulated_radiance);
						imgInterpolated.addColor(x, y, accumulated_radiance2);
					}
				}
			}
		}

		std::cout << std::endl;
	}

private:
	/******************************************************************
	* Check for closest intersection of a ray with the scene;
	* Returns true if intersection is found, as well as ray parameter
	* of intersection and id of intersected object
	*******************************************************************/
	bool Intersect_Scene(const Ray &ray, double *t, int *id, Vector *normal)
	{
		const int n = mTriangles.size();
		*t = 1e20;
		*id = -1;

		for (int i = 0; i < n; i++)
		{
			double d;
			if (!mTriangles[i].intersects(ray, d))
				continue;

			if (d > 0.0 && d < *t)
			{
				*t = d;
				*id = i;
				*normal = mTriangles[i].getNormal();
			}
		}
		return *t < 1e20;
	}

	/******************************************************************
	* Determine all form factors for all pairs of patches (of all
	* rectangles);
	* Evaluation of integrals in form factor equation is done via
	* Monte Carlo integration; mSamples are uniformly distributed and
	* equally weighted;
	* Computation accelerated by exploiting symmetries of form factor
	* estimation;
	*******************************************************************/
	void Calculate_Form_Factors(const int divisions, const int mc_sample)
	{
		/* Total number of patches in scene */
		const auto n = mTriangles.size();
		for (auto i = 0u; i < n; i++)
		{
			mTriangles[i].init_patchs(divisions);
			mPatchCount += mTriangles[i].getSubTriangleCount();
		}

		std::cout << "Number of triangles: " << n << std::endl;
		std::cout << "Number of patches: " << mPatchCount << std::endl;
		int mFormFactor_num = mPatchCount * mPatchCount;
		std::cout << "Number of form factors: " << mFormFactor_num << std::endl;

		/* 1D-array to hold form factor pairs */
		mFormFactor.clear();
		mFormFactor.resize(mFormFactor_num, 0.0);

		/* 1D-array with patch areas */
		std::vector<double> patch_area(mPatchCount, 0.0);

		/* Offsets for indexing of patches in 1D-array */
		std::vector<int> offset(n);

		for (auto i = 0; i < n; i++)
		{
			offset[i] = 0;
			for (auto k = 0; k < i; k++)
				offset[i] += mTriangles[k].getSubTriangleCount();
		}

		/* Precompute patch areas, assuming same size for each rectangle */
		for (auto i = 0; i < n; i++)
		{
			auto patch_i = offset[i];

			for (auto iSub = 0; iSub < mTriangles[i].getSubTriangleCount(); iSub++)
			{
				patch_area[patch_i + iSub] = mTriangles[i].getSubTriangle(iSub).area();
			}
		}


		/* Loop over all rectangles in scene */
		for (auto i = 0; i < n; i++)
		{
			auto patch_i = offset[i];

			std::cout << i << " ";

			/* Loop over all patches in rectangle i */
			for (auto iSub = 0; iSub < mTriangles[i].getSubTriangleCount(); iSub++)
			{
				std::cout << "*" << std::flush;

				const auto normal_i = mTriangles[i].getNormal();

				auto patch_j = 0;

				/* Loop over all rectangles in scene for rectangle i */
				for (auto j = 0; j < n; j++)
				{
					const auto normal_j = mTriangles[j].getNormal();

					/* Loop over all patches in rectangle j */
					for (auto jSub = 0; jSub < mTriangles[j].getSubTriangleCount(); jSub++)
					{
						/* Do not compute form factors for patches on same rectangle;
						also exploit symmetry to reduce computation;
						intemediate values; will be divided by patch area below */
						if (i < j)
						{
							double F = 0;

							/* Monte Carlo integration of form factor double integral */
							const auto Ni = mc_sample, Nj = mc_sample;

							/* Uniform PDF for Monte Carlo (1/Ai)x(1/Aj) */
							const auto pdf =
								(1.0 / patch_area[offset[i] + iSub]) *
								(1.0 / patch_area[offset[j] + jSub]);

							for (auto ias = 0; ias < Ni; ias++)
							{
								for (auto ibs = 0; ibs < Ni; ibs++)
								{
									for (auto jas = 0; jas < Nj; jas++)
									{
										for (auto jbs = 0; jbs < Nj; jbs++)
										{
											const auto xi = mTriangles[i].getSubTriangle(iSub).point_inside();
											const auto xj = mTriangles[j].getSubTriangle(jSub).point_inside();

											/* Check for visibility between sample points */
											const auto ij = (xj - xi).Normalized();

											double t;
											int id;
											Vector normal;
											if (Intersect_Scene(Ray(xi, ij), &t, &id, &normal) &&
												id != j)
											{
												continue; /* If intersection with other rectangle */
											}

											/* Cosines of angles beteen normals and ray inbetween */
											const auto d0 = normal_i.Dot(ij);
											const auto d1 = normal_j.Dot(-1.0 * ij);

											/* Continue if patches facing each other */
											if (d0 > 0.0 && d1 > 0.0)
											{
												/* Sample form factor */
												const auto K = d0 * d1 /
													(M_PI * (xj - xi).LengthSquared());

												/* Add weighted sample to estimate */
												F += K / pdf;
											}
										}
									}
								}
							}

							/* Divide by number of mSamples */
							F /= Ni* Ni* Nj* Nj;

							mFormFactor[patch_i * mPatchCount + patch_j] = F;
						}	
						
						patch_j++;			
					}
				}
				
				patch_i++;
			}

			std::cout << std::endl;
		}

		/* Copy upper to lower triangular values */
		for (auto i = 0u; i < mPatchCount - 1; i++)
		{
			for (auto j = i + 1; j < mPatchCount; j++)
			{
				mFormFactor[j * mPatchCount + i] = mFormFactor[i * mPatchCount + j];
			}
		}

		/* Divide by area to get final form factors */
		for (auto i = 0u; i < mPatchCount; i++)
		{
			for (auto j = 0u; j < mPatchCount; j++)
			{
				mFormFactor[i * mPatchCount + j] /= patch_area[i];

				/* Clamp to [0,1] */
				if (mFormFactor[i * mPatchCount + j] > 1.0)
					mFormFactor[i * mPatchCount + j] = 1.0;
			}
		}
	}

	/******************************************************************
	* Iterative computation of radiosity via Gathering; i.e. solution
	* using Gauss-Seidel iteration - reuse already computed values;
	* run-time O(n^2)
	*******************************************************************/
	void Calculate_Radiosity(const int iteration)
	{
		const int n = mTriangles.size();
		auto patch_i = 0;

		for (auto i = 0u; i < n; i++)
		{
			for (auto iSub = 0u; iSub < mTriangles[i].getSubTriangleCount(); iSub++)
			{
				Color B;

				auto patch_j = 0;
				for (auto j = 0u; j < n; j++)
				{
					for (auto jSub = 0u; jSub < mTriangles[j].getSubTriangleCount(); jSub++)
					{
						const auto Fij = mFormFactor[patch_i * mPatchCount + patch_j];

						/* Add form factor multiplied with radiosity of previous step */
						if (Fij > 0.0)
							B = B + Fij * mTriangles[j].getSubTriangle(jSub).getColor();

						patch_j++;
					}
				}


				/* Multiply sum with color of patch and add emission */
				B = mTriangles[i].getColor().MultComponents(B) + mTriangles[i].getEmission();

				/* Store overall patch radiosity of current iteration */
				mTriangles[i].getSubTriangle(iSub).setColor(B);

				patch_i++;
			}
		}
	}

	/******************************************************************
	* Compute radiance from radiosity by shooting rays into the scene;
	* Radiance directly proportional to radiosity for assumed diffuse
	* emitters/surfaces (multiply by PI);
	* At intersections either constant patch color is returned or a
	* smoothly interpolated color of 4x4 neighboring patches
	*******************************************************************/
	Color Radiance(unsigned int id, double distance, const Ray& ray, const int depth, bool interpolation = true)
	{
		const auto& tri = mTriangles[id];

		auto idSub = 0u;
		if (!tri.intersectsSub(ray, distance, idSub))
			assert(false && "No intersection");

		const auto& subTri = tri.getSubTriangle(idSub);

		/* Bicubic interpolation for smooth image */
		if (interpolation)
		{
			Triangle neigbours[4];
			tri.getNeighbours(idSub, neigbours);

			Color c[4];
			c[0] = neigbours[0].getColor();
			c[1] = neigbours[1].getColor();
			c[2] = neigbours[2].getColor();
			c[3] = neigbours[3].getColor();

			return cubicInterpolate(c, 0.5f) * Over_M_PI;
		}
		else
		{
			return subTri.getColor() * Over_M_PI;
		}
	}

};

#endif