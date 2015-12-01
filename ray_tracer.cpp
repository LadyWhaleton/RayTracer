/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;

    // TODO: determine the color
    // http://www.uio.no/studier/emner/matnat/ifi/INF3320/h03/undervisningsmateriale/lecture5.pdf
    // https://steveharveynz.wordpress.com/category/programming/c-raytracer/
    
    // I = I_ambient + I_diffuse + I_specular
    int numLights = world.lights.size();
    
    for (int i = 0; i < numLights; ++i)
    {
		Light* currLight = world.lights[i];
		Vector_3D<double> currLightColor = currLight->Emitted_Light(ray);
		
		// ======================================================
		// compute ambient component
		// I_a = ambientCoefficient * lightColor
		// ======================================================
		Vector_3D <double> ambient = color_ambient * currLightColor;
		
		// ======================================================
		// compute the diffuse component
		// I_d = diffuseCoefficient * lightColor * max (0, l.n)
		// ======================================================
		
		// l is the distance from the light source to the surface
		Vector_3D<double> l = currLight->position - intersection_point;
		l.Normalize();
		
		// I guess same_side_normal is surface normal, N
		double lDotN = Vector_3D<double>::Dot_Product(l, same_side_normal); 
		Vector_3D <double> diffuse = color_diffuse * currLightColor * max (0.0, lDotN);
		
		// ======================================================
		// compute the specular component.
		// I_s = specCoefficient * lightColor * max(0, v.r)^s
		// http://ogldev.atspace.co.uk/www/tutorial19/tutorial19.html
		// ======================================================
		
		// s is specular_power from Phong Shader class
		// v is distance from ray origin/endpoint to surface point, the surface to camera
		Vector_3D <double> v = ray.endpoint - intersection_point;
		v.Normalize();
		
		// http://www.gameprogrammer.net/delphi3dArchive/phongfordummies.htm
		// http://140.129.20.249/~jmchen/cg/docs/rendering%20pipeline/rendering/light_specular.html
		// R = 2(N.L)N - L
		Vector_3D <double> r = (same_side_normal * lDotN)*2 - l ; 
		r.Normalize();
		
		double vDotR = Vector_3D<double>::Dot_Product(v, r);
		double maxVal = max(0.0, vDotR);
		Vector_3D <double> specular = color_specular * currLightColor * pow( maxVal, specular_power);
		
		color += diffuse + specular;
		
	}

    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;
    // TODO: determine the color

    return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
    // TODO
    // http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    // https://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter1.htm
    // http://www.sfdm.scad.edu/faculty/mkesson/vsfx419/wip/best/winter12/jonathan_mann/raytracer.html
    // a = D^2, b = 2D(O-C), c = |O - C|^2 - r^2
    double a = Vector_3D <double>::Dot_Product(ray.direction, ray.direction);
    double b = 2 * Vector_3D <double>:: Dot_Product(ray.direction, ray.endpoint - center); 
    double c = Vector_3D <double>::Dot_Product(ray.endpoint - center, ray.endpoint - center) - radius*radius; 
    double discr = b*b - 4 * a * c;

    if (discr < 0)     // if discriminant < 0, no roots
        return false;
    
    else
    {
		double t1, t2, t; // these represent the solutions for t
		t1 = -0.5 * (b + sqrt(discr)) / a;
		t2 = -0.5 * (b - sqrt(discr)) / a;
		
		// find the smaller positive root and set that to t_max
		if (t1 >= 0 && t2 < 0) // t1 is pos
			t = t1;
		else if (t1 < 0 && t2 >= 0) // t2 is pos
			t = t2;
		else if (t1 >= 0 && t2 >= 0)
			t = min(t1, t2);
        else 
            return false;
		
        // t has to be bigger than small_t to register an intersection with a ray	
		if (t > small_t)
        {
			ray.semi_infinite = false;
            ray.t_max = t;
            ray.current_object = this;  // set current object
            return true;
        }		
	}

    return false;       
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal = location - center;
    normal.Normalize();
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
    // TODO
    // https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld017.htm
    // http://www.prinmath.com/csci5229/Sp11/handouts/obj-ray-inter.pdf
    // (a - o) * n / (d * n), o = endpoint, d = direction, n = normal
    // x1 = a, which is some point on the plane
    double top = Vector_3D<double>::Dot_Product(x1 - ray.endpoint, normal);
    double bot = Vector_3D<double>::Dot_Product(ray.direction, normal);
    double t = top / bot;
    
    // if t > 0, ray towards plane and will eventually intersect it
    // t has to be bigger than small_t to register an intersection with a ray
    if (t > small_t)
    {
		ray.t_max = t;
		ray.current_object = this;
		ray.semi_infinite = false;
		return true;
	}
    
	// if t < 0, ray is away from plane and will never intersect	
    return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
    // TODO 
    // We are essentially converting screen/camera space to world space
 
    // Need to get the grid location of pixel_index (x and y coords)
    Vector_2D<double> gridPos = film.pixel_grid.X(pixel_index);
    
    // convert camera/screen coordinates (gridPosition) to world coordinates
    Vector_3D<double> worldPos = horizontal_vector*gridPos.x + vertical_vector*gridPos.y;
    
    // focal point = eye vector (the origin) / the image plane
    Vector_3D<double> result = focal_point + worldPos;
    
    return result;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    //cout << "Looking for closest intersection!\n";

    // TODO: start
    // create a variable to hold the current closest object (initially infinity)
    int numObjects = objects.size();

    // for each object in the scene
    for (int i=0; i < numObjects; i++)
        objects[i]->Intersection(ray);

    return ray.current_object;
}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO: start
    // set up the initial view ray here
    Ray ray (camera.position, camera.World_Position(pixel_index) - camera.position);
	ray.t_max = FLT_MAX;
    ray.semi_infinite = true;
	
    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
    // TODO: start
    Vector_3D<double> color;
    const Object* closestObj = Closest_Intersection(ray);

    if (closestObj != NULL)
    {
        Vector_3D <double> pt = ray.Point(ray.t_max);
        Vector_3D <double> n = closestObj->Normal(pt);
        return closestObj->material_shader->Shade_Surface(ray, *closestObj, pt, n);
    }

    return color;
}
