#ifndef RendererRectangles_H_included
#define RendererRectangles_H_included

#include <vector>

#include "Renderer.hpp"
#include "Rectangle.hpp"

#define drand48() (((double)rand())/((double)RAND_MAX))

class RendererRectangles : public Renderer
{
	const double Over_M_PI = 1.0 / M_PI;

	std::vector<double> mFormFactor;
	unsigned int mPatchCount;
	std::vector<Rectangle> mRectangles;
	

public:
	RendererRectangles(unsigned int mWidth, unsigned int mHeight, unsigned int mSamples) : Renderer(mWidth, mHeight, mSamples)
	{
		mPatchCount = 0;
	}

	virtual void buildScene(const std::vector<Rectangle>& rects)
	{
		mRectangles = rects;
	}

	virtual void Render(Image& img, Image& imgInterpolated)
	{
		/* Set camera origin and viewing direction (negative z direction) */
		Ray camera(Vector(50.0, 52.0, 295.6), Vector(0.0, -0.042612, -1.0).Normalized());

		/* Image edge vectors for pixel sampling */
		Vector cx = Vector(mWidth * 0.5135 / mHeight);
		Vector cy = (cx.Cross(camera.dir)).Normalized() * 0.5135;

		std::cout << "Calculating form factors" << std::endl;
		int patches_a = 4;
		int patches_b = 4;
		int MC_mSamples = 3;

		Calculate_Form_Factors(patches_a, patches_b, MC_mSamples);

		/* Iterative solution of radiosity linear system */
		std::cout << "Calculating radiosity" << std::endl;
		int iterations = 40;
		for (int i = 0; i < iterations; i++)
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
						Color accumulated_radiance = Color();
						Color accumulated_radiance2 = Color();

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

							/* Determine constant radiance */
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
			const int n = mRectangles.size();
			*t = 1e20;
			*id = -1;

			for (int i = 0; i < n; i++)
			{
				double d = mRectangles[i].intersect(ray);
				if (d > 0.0 && d < *t)
				{
					*t = d;
					*id = i;
					*normal = mRectangles[i].normal;
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
		void Calculate_Form_Factors(const int a_div_num, const int b_div_num, const int mc_sample)
		{
			/* Total number of patches in scene */
			const int n = mRectangles.size();
			for (int i = 0; i < n; i++)
			{
				mRectangles[i].init_patchs(a_div_num, b_div_num);
				mPatchCount += mRectangles[i].a_num * mRectangles[i].b_num;
			}

			std::cout << "Number of rectangles: " << n << std::endl;
			std::cout << "Number of patches: " << mPatchCount << std::endl;
			int mFormFactor_num = mPatchCount * mPatchCount;
			std::cout << "Number of form factors: " << mFormFactor_num << std::endl;

			/* 1D-array to hold form factor pairs */
			mFormFactor.clear();
			mFormFactor.resize(mFormFactor_num);
			memset(mFormFactor.data(), 0, sizeof(double)* mFormFactor_num);

			/* 1D-array with patch areas */
			double *patch_area = new double[mPatchCount];
			memset(patch_area, 0, sizeof(double)* mPatchCount);

			/* Precompute patch areas, assuming same size for each rectangle */
			for (int i = 0; i < n; i++)
			{
				int patch_i = 0;
				for (int k = 0; k < i; k++)
					patch_i += mRectangles[k].a_num * mRectangles[k].b_num;

				for (int ia = 0; ia < mRectangles[i].a_num; ia++)
				{
					for (int ib = 0; ib < mRectangles[i].b_num; ib++)
					{
						patch_area[patch_i + ia* mRectangles[i].b_num + ib] =
							((mRectangles[i].edge_a / mRectangles[i].a_num).
							Cross((mRectangles[i].edge_b / mRectangles[i].b_num))).Length();
					}
				}
			}

			/* Offsets for indexing of patches in 1D-array */
			int *offset = new int[n];

			for (int i = 0; i < n; i++)
			{
				offset[i] = 0;
				for (int k = 0; k < i; k++)
					offset[i] += mRectangles[k].a_num * mRectangles[k].b_num;
			}

			/* Loop over all rectangles in scene */
			for (int i = 0; i < n; i++)
			{
				int patch_i = offset[i];

				std::cout << i << " ";

				/* Loop over all patches in rectangle i */
				for (int ia = 0; ia < mRectangles[i].a_num; ia++)
				{
					std::cout << "*" << std::flush;
					for (int ib = 0; ib < mRectangles[i].b_num; ib++)
					{
						const Vector normal_i = mRectangles[i].normal;

						int patch_j = 0;

						/* Loop over all rectangles in scene for rectangle i */
						for (int j = 0; j < n; j++)
						{
							const Vector normal_j = mRectangles[j].normal;

							/* Loop over all patches in rectangle j */
							for (int ja = 0; ja < mRectangles[j].a_num; ja++)
							{
								for (int jb = 0; jb < mRectangles[j].b_num; jb++)
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
											(1.0 / patch_area[offset[i] + ia*mRectangles[i].b_num + ib]) *
											(1.0 / patch_area[offset[j] + ja*mRectangles[j].b_num + jb]);

										/* Determine rays of NixNi uniform mSamples of patch
										on i to NjxNj uniform mSamples of patch on j */
										for (int ias = 0; ias < Ni; ias++)
										{
											for (int ibs = 0; ibs < Ni; ibs++)
											{
												for (int jas = 0; jas < Nj; jas++)
												{
													for (int jbs = 0; jbs < Nj; jbs++)
													{
														/* Determine sample points xi, xj on both patches */
														const double u0 = (double)(ias + 0.5) / Ni,
															u1 = (double)(ibs + 0.5) / Ni;
														const double u2 = (double)(jas + 0.5) / Nj,
															u3 = (double)(jbs + 0.5) / Nj;

														const Vector xi = mRectangles[i].p0 +
															mRectangles[i].edge_a * ((double)(ia + u0) / mRectangles[i].a_num) +
															mRectangles[i].edge_b * ((double)(ib + u1) / mRectangles[i].b_num);
														const Vector xj = mRectangles[j].p0 +
															mRectangles[j].edge_a * ((double)(ja + u2) / mRectangles[j].a_num) +
															mRectangles[j].edge_b * ((double)(jb + u3) / mRectangles[j].b_num);

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

										/* Divide by number of mSamples */
										F /= (Ni)* (Ni)* (Nj)* (Nj);

										mFormFactor[patch_i * mPatchCount + patch_j] = F;
									}
									patch_j++;
								}
							}
						}
						patch_i++;
					}
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
			const int n = mRectangles.size();
			int patch_i = 0;

			for (int i = 0; i < n; i++)
			{
				for (int ia = 0; ia < mRectangles[i].a_num; ia++)
				{
					for (int ib = 0; ib < mRectangles[i].b_num; ib++)
					{
						Color B;

						int patch_j = 0;
						for (int j = 0; j < n; j++)
						{
							for (int ja = 0; ja < mRectangles[j].a_num; ja++)
							{
								for (int jb = 0; jb < mRectangles[j].b_num; jb++)
								{
									const double Fij = mFormFactor[patch_i * mPatchCount + patch_j];

									/* Add form factor multiplied with radiosity of previous step */
									if (Fij > 0.0)
										B = B + Fij * mRectangles[j].patch[ja * mRectangles[j].b_num + jb];

									patch_j++;
								}
							}
						}
						/* Multiply sum with color of patch and add emission */
						B = mRectangles[i].color.MultComponents(B) + mRectangles[i].emission;

						/* Store overall patch radiosity of current iteration */
						mRectangles[i].patch[ia * mRectangles[i].b_num + ib] = B;
						patch_i++;
					}
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
		Color Radiance(int id, double t, Ray ray, const int depth, bool interpolation = true)
		{
			/* Determine intersection point on rectangle */
			const Rectangle &obj = mRectangles[id];
			const Vector hitpoint = ray.org + t * ray.dir;

			/* Determine intersected patch */
			const Vector v = hitpoint - obj.p0;
			const double a_len = v.Dot(obj.edge_a.Normalized());
			const double b_len = v.Dot(obj.edge_b.Normalized());

			double da = obj.a_num * a_len / obj.a_len;
			double db = obj.b_num * b_len / obj.b_len;

			int ia = int(da); if (ia >= obj.a_num) ia--;
			int ib = int(db); if (ib >= obj.b_num) ib--;

			/* Bicubic interpolation for smooth image */
			if (interpolation)
			{
				Color c[4][4];

				int ia = int(da - 0.5);
				int ib = int(db - 0.5);

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
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

};

#endif