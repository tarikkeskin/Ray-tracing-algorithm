#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include<float.h>


using namespace parser;
using namespace std;

typedef unsigned char RGB[3];

struct Ray
{
	Vec3f origin;
	Vec3f direction;
	Vec3f intersectionPoint;
	Vec3f normalVector;
	bool shadowRay;
	bool hitFlag;
	int material_id;
	float t;
};

Vec3f subtract(const Vec3f &a, const Vec3f &b)
{
	Vec3f temp;
	temp.x = a.x - b.x;
	temp.y = a.y - b.y;
	temp.z = a.z - b.z;
	return temp;
}


Vec3f add(const Vec3f &a, const Vec3f &b)
{
	Vec3f temp;
	temp.x = a.x + b.x;
	temp.y = a.y + b.y;
	temp.z = a.z + b.z;
	return temp;
}

Vec3f mul(const Vec3f &a, const Vec3f &b)
{
	Vec3f temp;
	temp.x = a.x * b.x;
	temp.y = a.y * b.y;
	temp.z = a.z * b.z;
	return temp;
}

Vec3f mulcons(const Vec3f &a, float b)
{
	Vec3f temp;
	temp.x = a.x * b;
	temp.y = a.y * b;
	temp.z = a.z * b;
	return temp;
}


float dotProduct(const Vec3f &a, const Vec3f &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}


Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.y*b.z - a.z*b.y;
	result.y = a.z*b.x - a.x*b.z;
	result.z = a.x*b.y - a.y*b.x;

	return result;
}


float findLength(const Vec3f &a)
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}


float findDistance(const Vec3f &a, const Vec3f &b)
{
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2)
{
	return v0.x*(v1.y*v2.z - v2.y*v1.z) + v0.y*(v2.x*v1.z - v1.x*v2.z) + v0.z * (v1.x*v2.y - v1.y*v2.x);
}


Vec3f normalizeVector(const Vec3f &a)
{
	Vec3f result;
	result.x = a.x / findLength(a);
	result.y = a.y / findLength(a);
	result.z = a.z / findLength(a);
	return result;
}


Vec3f intersectPoint(const Ray &ray, float t)
{
	Vec3f result;
	result.x = ray.origin.x + t * ray.direction.x;
	result.y = ray.origin.y + t * ray.direction.y;
	result.z = ray.origin.z + t * ray.direction.z;

	return result;
}

Ray generateRay(const Camera &camera, int i, int j)
{
	Vec3f m, q, u, v, s;

	float su = (camera.near_plane.y - camera.near_plane.x)*(j + 0.5) / camera.image_width;
	float sv = (camera.near_plane.w - camera.near_plane.z)*(i + 0.5) / camera.image_height;

	//m = e + -wdistance
	m.x = camera.position.x + (camera.gaze.x * camera.near_distance);
	m.y = camera.position.y + (camera.gaze.y * camera.near_distance);
	m.z = camera.position.z + (camera.gaze.z * camera.near_distance);

	u = crossProduct(camera.gaze, camera.up);
	u = normalizeVector(u);

	v = crossProduct(u, camera.gaze);
	v = normalizeVector(v);

	//q = m + lu + tv
	q.x = m.x + (u.x*camera.near_plane.x) + (v.x*camera.near_plane.w);
	q.y = m.y + (u.y*camera.near_plane.x) + (v.y*camera.near_plane.w);
	q.z = m.z + (u.z*camera.near_plane.x) + (v.z*camera.near_plane.w);

	//s = q + su.u - sv.v
	s.x = q.x + (u.x*su) - (v.x * sv);
	s.y = q.y + (u.y*su) - (v.y * sv);
	s.z = q.z + (u.z*su) - (v.z * sv);

	Vec3f origin = camera.position;
	Vec3f direction = normalizeVector(subtract(s, camera.position));
	//r(t) = e + (s – e)t = e + dt
	Ray ray{ origin,direction,false };

	return ray;
}


Ray sphereIntersection(const Ray &ray, const Scene &scene, const Sphere &sphere)
{
	Ray hit;

	const float a = dotProduct(ray.direction, ray.direction);
	Vec3f v1 = subtract(ray.origin, scene.vertex_data[sphere.center_vertex_id - 1]);
	const float b = 2 * dotProduct(ray.direction, v1);
	const float c = dotProduct(v1, v1) - (sphere.radius * sphere.radius);

	if (b*b - 4 * a*c < 0)	hit.hitFlag = false;
	else
	{
		const float t1 = (-1 * b + sqrtf(b*b - 4 * a*c)) / 2 * a;
		const float t2 = (-1 * b - sqrtf(b*b - 4 * a*c)) / 2 * a;

		hit.material_id = sphere.material_id;
		hit.hitFlag = true;

		const float t = fmin(t1, t2);
		hit.intersectionPoint = intersectPoint(ray, t);
		hit.normalVector = subtract(hit.intersectionPoint, scene.vertex_data[sphere.center_vertex_id - 1]);
		hit.normalVector.x /= sphere.radius;
		hit.normalVector.y /= sphere.radius;
		hit.normalVector.z /= sphere.radius;

		hit.t = t;
	}

	return hit;
}

Ray triangleIntersection(const Ray &ray, const Vec3f &a, const Vec3f &b, const Vec3f &c, int material_id)
{
	Ray hit;
	hit.hitFlag = false;

	Vec3f o = ray.origin;
	Vec3f d = ray.direction;

	Vec3f v1 = subtract(a, b);
	Vec3f v2 = subtract(a, c);
	Vec3f v3 = subtract(a, o);

	float detA = determinant(v1, v2, d);
	if (detA == 0.0) return hit;

	float t = (determinant(v1, v2, v3)) / detA;
	if (t <= 0.0) {
		return hit;
	}

	float gamma = (determinant(v1, v3, d)) / detA;
	if (gamma < 0 || gamma > 1) {
		return hit;
	}

	float beta = (determinant(v3, v2, d)) / detA;
	if (beta < 0 || beta >(1 - gamma)) {
		return hit;
	}

	hit.normalVector = normalizeVector(crossProduct(subtract(b, a), subtract(c, a)));
	hit.intersectionPoint = intersectPoint(ray, t);
	hit.hitFlag = true;
	hit.material_id = material_id;
	hit.t = t;



	return hit;
}


Ray meshIntersection(const Ray &ray, const Mesh &mesh, const Scene &scene, int material_id)
{
	Ray hit;
	hit.hitFlag = false;
	vector<Ray> rayInfo;

	for (int i = 0; i < mesh.faces.size(); i++)
	{
		Vec3f v0 = scene.vertex_data[mesh.faces[i].v0_id - 1];
		Vec3f v1 = scene.vertex_data[mesh.faces[i].v1_id - 1];
		Vec3f v2 = scene.vertex_data[mesh.faces[i].v2_id - 1];

		hit = triangleIntersection(ray, v0, v1, v2, material_id);
		if (hit.hitFlag && hit.t >= 0)
		{
			hit.normalVector = normalizeVector(crossProduct(subtract(v1, v0), subtract(v2, v0)));
			hit.intersectionPoint = intersectPoint(ray, hit.t);
			hit.material_id = material_id;

			rayInfo.push_back(hit);
		}
	}

	hit.hitFlag = false;

	if (rayInfo.size() != 0)
	{
		hit = rayInfo[0];

		for (int i = 1; i < rayInfo.size(); i++)
		{
			if (rayInfo[i].t < hit.t)
			{
				hit = rayInfo[i];
			}
		}
		hit.hitFlag = true;
	}
	return hit;
}


vector<Vec3f> BoundingVolumeCalculate(const Mesh mesh, const vector<Vec3f> vertex_data) {
	vector<Vec3f> vol;
	Vec3f v0 = vertex_data[mesh.faces[0].v0_id - 1];//initialize the min max coordinates (Can't initialize 0 0 0 because it may be outside of Volume)
	Vec3f v1, v2;
	Vec3f vec;

	float x_min = v0.x;
	float y_min = v0.y;
	float z_min = v0.z;
	float x_max = v0.x;
	float y_max = v0.y;
	float z_max = v0.z;

	for (int faceNumber = 0; faceNumber < mesh.faces.size(); faceNumber++)
	{
		v0 = vertex_data[mesh.faces[faceNumber].v0_id - 1];
		v1 = vertex_data[mesh.faces[faceNumber].v1_id - 1];
		v2 = vertex_data[mesh.faces[faceNumber].v2_id - 1];

		if (v0.x < x_min) { x_min = v0.x; }
		if (v0.y < y_min) { y_min = v0.y; }
		if (v0.z < z_min) { z_min = v0.z; }
		if (v0.x > x_max) { x_max = v0.x; }
		if (v0.y > y_max) { y_max = v0.y; }
		if (v0.z > z_max) { z_max = v0.z; }

		if (v1.x < x_min) { x_min = v1.x; }
		if (v1.y < y_min) { y_min = v1.y; }
		if (v1.z < z_min) { z_min = v1.z; }
		if (v1.x > x_max) { x_max = v1.x; }
		if (v1.y > y_max) { y_max = v1.y; }
		if (v1.z > z_max) { z_max = v1.z; }

		if (v2.x < x_min) { x_min = v2.x; }
		if (v2.y < y_min) { y_min = v2.y; }
		if (v2.z < z_min) { z_min = v2.z; }
		if (v2.x > x_max) { x_max = v2.x; }
		if (v2.y > y_max) { y_max = v2.y; }
		if (v2.z > z_max) { z_max = v2.z; }
	}

	vec.x = x_min;
	vec.y = y_min;
	vec.z = z_min;
	vol.push_back(vec);

	vec.x = x_max;
	vec.y = y_min;
	vec.z = z_min;
	vol.push_back(vec);

	vec.x = x_min;
	vec.y = y_min;
	vec.z = z_max;
	vol.push_back(vec);

	vec.x = x_max;
	vec.y = y_min;
	vec.z = z_max;
	vol.push_back(vec);

	vec.x = x_min;
	vec.y = y_max;
	vec.z = z_min;
	vol.push_back(vec);

	vec.x = x_max;
	vec.y = y_max;
	vec.z = z_min;
	vol.push_back(vec);

	vec.x = x_min;
	vec.y = y_max;
	vec.z = z_max;
	vol.push_back(vec);

	vec.x = x_max;
	vec.y = y_max;
	vec.z = z_max;
	vol.push_back(vec);

	return vol;

}

bool BoundingVolumeTest(Ray r, vector<Vec3f> vol) {
	bool hit_happened = false;
	if (triangleIntersection(r, vol[0], vol[1], vol[4], 0).hitFlag) { hit_happened = true; } //front face
	if (triangleIntersection(r, vol[5], vol[4], vol[1], 0).hitFlag) { hit_happened = true; }
	if (triangleIntersection(r, vol[0], vol[2], vol[1], 0).hitFlag) { hit_happened = true; } //down face
	if (triangleIntersection(r, vol[2], vol[3], vol[1], 0).hitFlag) { hit_happened = true; }
	if (triangleIntersection(r, vol[6], vol[7], vol[3], 0).hitFlag) { hit_happened = true; } // back face
	if (triangleIntersection(r, vol[6], vol[3], vol[2], 0).hitFlag) { hit_happened = true; }
	if (triangleIntersection(r, vol[7], vol[6], vol[5], 0).hitFlag) { hit_happened = true; } // up face
	if (triangleIntersection(r, vol[6], vol[4], vol[5], 0).hitFlag) { hit_happened = true; }
	if (triangleIntersection(r, vol[4], vol[2], vol[0], 0).hitFlag) { hit_happened = true; } // left face
	if (triangleIntersection(r, vol[6], vol[2], vol[4], 0).hitFlag) { hit_happened = true; }
	if (triangleIntersection(r, vol[7], vol[5], vol[3], 0).hitFlag) { hit_happened = true; } //right face
	if (triangleIntersection(r, vol[5], vol[1], vol[3], 0).hitFlag) { hit_happened = true; }
	return hit_happened;
}

Vec3f Irradiance(const PointLight &light, const Vec3f &intersectionPoint)
{
	Vec3f irradiance;
	Vec3f v = subtract(light.position, intersectionPoint);
	float vq = dotProduct(v, v);

	if (vq == 0.0)
	{

		irradiance.x = 0;
		irradiance.y = 0;
		irradiance.z = 0;
	}
	else {
		irradiance.x = light.intensity.x / vq;
		irradiance.y = light.intensity.y / vq;
		irradiance.z = light.intensity.z / vq;
	}
	return irradiance;
}


const Vec3f Diffuse(const PointLight &light, const Scene &scene, int material_id, const Vec3f &normal, const Vec3f &intersectionPoint)
{
	Vec3f diffuse;

	Vec3f irradiance = Irradiance(light, intersectionPoint);

	Vec3f l = subtract(light.position, intersectionPoint);
	l = normalizeVector(l);

	float dot = dotProduct(l, normal);
	if (dot < 0) dot = 0;

	Vec3f temp;
	temp = mulcons(scene.materials[material_id - 1].diffuse, dot);
	diffuse = mul(temp, irradiance);
	return diffuse;
}


Vec3f Specular(const PointLight &light, const Scene &scene, const Ray &ray, int material_id, const Vec3f &normal, const Vec3f &intersectionPoint)
{
	Vec3f specular;

	Material material = scene.materials[material_id - 1];

	Vec3f irradiance = Irradiance(light, intersectionPoint);

	Vec3f w = subtract(light.position, intersectionPoint);
	w = normalizeVector(w);

	Vec3f n = subtract(w, ray.direction);
	n = normalizeVector(n);

	float dot = dotProduct(normal, n);
	if (dot < 0) dot = 0;

	specular = mul(material.specular, irradiance);
	specular = mulcons(specular, pow(dot, material.phong_exponent));

	return specular;
}

bool shadowEachObject(const Scene &scene, const Ray &shadowRay, bool shadowFlag, float x_distance, vector<vector<Vec3f>> bounding_mesh) {

	Ray ray;
	for (int sphereNumber = 0; sphereNumber < scene.spheres.size(); sphereNumber++)
	{
		ray = sphereIntersection(shadowRay, scene, scene.spheres[sphereNumber]);

		if (ray.hitFlag)
		{
			if (x_distance>ray.t && ray.t>=0)
			{
				shadowFlag = true;
			}
		}
	}
	//triangles
	for (int triangleNumber = 0; triangleNumber < scene.triangles.size(); triangleNumber++)
	{
		Triangle currentTriangle = scene.triangles[triangleNumber];
		Vec3f v0 = scene.vertex_data[currentTriangle.indices.v0_id - 1];
		Vec3f v1 = scene.vertex_data[currentTriangle.indices.v1_id - 1];
		Vec3f v2 = scene.vertex_data[currentTriangle.indices.v2_id - 1];

		ray = triangleIntersection(shadowRay, v0, v1, v2, currentTriangle.material_id);

		if (ray.hitFlag)
		{
			if (x_distance>ray.t && ray.t>=0)
			{
				shadowFlag = true;
			}
		}
	}

	//meshes
	for (int meshNumber = 0; meshNumber < scene.meshes.size(); meshNumber++)
	{
		Mesh currentMesh = scene.meshes[meshNumber];

		if (BoundingVolumeTest(shadowRay, bounding_mesh[meshNumber])) {
			ray = meshIntersection(shadowRay, currentMesh, scene, currentMesh.material_id);
		}
		if (ray.hitFlag)
		{
			if (x_distance>ray.t && ray.t>=0)
			{
				shadowFlag = true;
			}
		}
	}
	return shadowFlag;
}



Vec3f shading(const Scene &scene, const Ray &res, const Camera &currentCamera, const Ray &ray, int maxDepth ,const vector<vector<Vec3f>> bounding_mesh)
{
	Vec3f pixel, color;

	pixel.x = scene.background_color.x;
	pixel.y = scene.background_color.y;
	pixel.z = scene.background_color.z;

	int material_id = res.material_id;

	pixel = mul(scene.materials[material_id - 1].ambient, scene.ambient_light);


	for (int lightNumber = 0; lightNumber < scene.point_lights.size(); lightNumber++)
	{
		bool shadowFlag = false;

		PointLight currentLight = scene.point_lights[lightNumber];
		float lightdirection = findDistance(currentLight.position, currentCamera.position);

		Vec3f wi = subtract(currentLight.position, res.intersectionPoint);
		wi = normalizeVector(wi);

		Vec3f epsilonvector;
		epsilonvector = mulcons(wi, scene.shadow_ray_epsilon);

		Vec3f orig;
		orig = add(res.intersectionPoint, epsilonvector);

		Ray shadowRay = { orig,wi,true };

		float x_distance = subtract(currentLight.position, shadowRay.origin).x / shadowRay.direction.x;

		shadowFlag = shadowEachObject(scene, shadowRay, shadowFlag, x_distance, bounding_mesh);

		if (!shadowFlag || lightdirection == 0.0)
		{
			int material_id = res.material_id;

			Vec3f diffuse = Diffuse(currentLight, scene, material_id, res.normalVector, res.intersectionPoint);

			Vec3f specular = Specular(currentLight, scene, ray, material_id, res.normalVector, res.intersectionPoint);

			Vec3f difspecular;
			difspecular = add(specular, diffuse);
			pixel = add(pixel, difspecular);
		}

	}
	Vec3f reflection;

	reflection.x = 0;
	reflection.y = 0;
	reflection.z = 0;

	if (maxDepth > 0 && (scene.materials[material_id - 1].mirror.x > 0 || scene.materials[material_id - 1].mirror.y > 0 || scene.materials[material_id - 1].mirror.z > 0))
	{
		float d = -2 * dotProduct(ray.direction, res.normalVector);

		Vec3f wi2, v;
		v = mulcons(res.normalVector, d);
		wi2 = add(v, ray.direction);

		wi2 = normalizeVector(wi2);

		Vec3f epsilon;

		epsilon = mulcons(wi2, scene.shadow_ray_epsilon);

		Ray reflectionRay = { add(res.intersectionPoint, epsilon),wi2,false };

		vector<Ray> rayInfo;

		for (int k = 0; k < scene.spheres.size(); k++) {

			Ray hit = sphereIntersection(reflectionRay, scene, scene.spheres[k]);

			if (hit.hitFlag && hit.t >= 0)
			{
				rayInfo.push_back(hit);
			}

		}
		for (int k = 0; k < scene.triangles.size(); k++) {

			Vec3f v0 = scene.vertex_data[scene.triangles[k].indices.v0_id - 1];
			Vec3f v1 = scene.vertex_data[scene.triangles[k].indices.v1_id - 1];
			Vec3f v2 = scene.vertex_data[scene.triangles[k].indices.v2_id - 1];

			Ray hit = triangleIntersection(reflectionRay, v0, v1, v2, scene.triangles[k].material_id);

			if (hit.hitFlag && hit.t >= 0)
			{
				rayInfo.push_back(hit);
			}
		}

		for (int k = 0; k < scene.meshes.size(); k++) {

			if (BoundingVolumeTest(reflectionRay, bounding_mesh[k])) {
				Ray hit = meshIntersection(reflectionRay, scene.meshes[k], scene, scene.meshes[k].material_id);

				if (hit.hitFlag && hit.t >= 0)
				{
					rayInfo.push_back(hit);
				}
			}
		}

		Ray res;
		res.hitFlag = false;

		if (rayInfo.size() != 0)
		{
			res = rayInfo[0];

			for (int i = 1; i < rayInfo.size(); i++)
			{
				if (rayInfo[i].t < res.t)
				{
					res = rayInfo[i];
				}
			}
			res.hitFlag = true;
		}

		if (res.hitFlag != false) {
			reflection = shading(scene, res, currentCamera, reflectionRay, maxDepth - 1, bounding_mesh);

			Vec3f temp;
			temp = mul(reflection, scene.materials[material_id - 1].mirror);
			pixel = add(pixel, temp);
		}
		else {
			Vec3f empty_pixel;
			empty_pixel.x = float(scene.background_color.x);
			empty_pixel.y = float(scene.background_color.y);
			empty_pixel.z = float(scene.background_color.z);
			pixel = add(pixel,empty_pixel);
		}
	}


	color = pixel;
	return color;
}




int main(int argc, char* argv[])
{
	// Sample usage for reading an XML scene file
	Scene scene;

	scene.loadFromXml(argv[1]);
	// bounding box
	vector<vector<Vec3f>> meshBoxes;

	for (int meshNumber = 0; meshNumber < scene.meshes.size(); meshNumber++) {
		meshBoxes.push_back( BoundingVolumeCalculate(scene.meshes[meshNumber], scene.vertex_data) );
	}

	// for each camera
	for (int cameraNumber = 0; cameraNumber < scene.cameras.size(); cameraNumber++)
	{
		Camera currentCamera = scene.cameras[cameraNumber];

		unsigned char* image = new unsigned char[currentCamera.image_width * currentCamera.image_height * 3];
		int imageN = 0;

		for (int i = 0; i < currentCamera.image_height; i++)
		{
			for (int j = 0; j < currentCamera.image_width; j++)
			{
				//cout << i << " " << j << endl;
				//for each pixel ray

				double tmin = DBL_MAX; // MAX DOUBLE
				int closestObj = -1;
				Ray ray = generateRay(currentCamera, i, j);

				vector<Ray> rayInfo;

				for (int k = 0; k < scene.spheres.size(); k++) {

					Ray hit = sphereIntersection(ray, scene, scene.spheres[k]);

					if (hit.hitFlag && hit.t >= 0)
					{
						rayInfo.push_back(hit);
					}

				}
				for (int k = 0; k < scene.triangles.size(); k++) {

					Vec3f v0 = scene.vertex_data[scene.triangles[k].indices.v0_id - 1];
					Vec3f v1 = scene.vertex_data[scene.triangles[k].indices.v1_id - 1];
					Vec3f v2 = scene.vertex_data[scene.triangles[k].indices.v2_id - 1];

					Ray hit = triangleIntersection(ray, v0, v1, v2, scene.triangles[k].material_id);

					if (hit.hitFlag && hit.t >= 0)
					{
						rayInfo.push_back(hit);
					}
				}

				for (int k = 0; k < scene.meshes.size(); k++) {

					if (BoundingVolumeTest(ray , meshBoxes[k])) {
						Ray hit = meshIntersection(ray, scene.meshes[k], scene, scene.meshes[k].material_id);

						if (hit.hitFlag && hit.t >= 0)
						{
							rayInfo.push_back(hit);

						}
					}
				}

				Ray res;
				res.hitFlag = false;

				if (rayInfo.size() != 0)
				{
					res = rayInfo[0];

					for (int i = 1; i < rayInfo.size(); i++)
					{
						if (rayInfo[i].t < res.t)
						{
							res = rayInfo[i];
						}
					}
					res.hitFlag = true;
				}

				if (res.hitFlag)
				{
					Vec3f color = shading(scene, res, currentCamera, ray, (scene.max_recursion_depth), meshBoxes);

					if (color.x < 0) image[imageN++] = 0;
					if (color.x > 255) image[imageN++] = 255;
					else image[imageN++] = round(color.x);

					if (color.y < 0) image[imageN++] = 0;
					if (color.y > 255) image[imageN++] = 255;
					else image[imageN++] = round(color.y);

					if (color.z < 0) image[imageN++] = 0;
					if (color.z > 255) image[imageN++] = 255;
					else image[imageN++] = round(color.z);

				}
				else
				{
					image[imageN++] = scene.background_color.x;
					image[imageN++] = scene.background_color.y;
					image[imageN++] = scene.background_color.z;
				}
			}
		}

		write_ppm(currentCamera.image_name.c_str(), image, currentCamera.image_width, currentCamera.image_height);
	}

}
