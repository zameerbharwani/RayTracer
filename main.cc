#include <iostream>
#include "color.h"
#include "vec3.h"
#include "ray.h"

using namespace std;

bool hit_sphere(const point3& center, double radius, const ray& r) {

    /*
       A := point of origin of the ray 
       b := direction vector of the ray
       C := center of the sphere

       0 = t^2 * b^2 + 2*t*b*(A-C) + (A-C)*(A-C) - r^2
       coeff_1  = b^2
       coeff_2 = 2*b*(A-C) 
       coeff_3 = (A-C)*(A-C) - r^2
       Quadratic in terms of t
    */
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius * radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

color ray_color(const ray &r) {
    
    if (hit_sphere(point3(0,0,-1),0.5,r)) return color(1,0,0);

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0,1.0,1.0) + t*color(0.5,0.7,1.0);
}

int main() {
    
    // Image
    const auto aspect_ratio = 16.0/9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Camera
    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = point3(0,0,0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0,0,focal_length);

    // Render
    cout<<"P3"<<endl<<image_width<< ' '<<image_height<<endl<<"255"<<endl;

    for (int j = image_height -1; j>=0; --j) {
        cerr<<"\rScanlines remaining: " << j << ' '<<flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / (image_width-1);
            auto v = double(j) / (image_height-1);
            ray r(origin, lower_left_corner + u*horizontal + v*vertical - origin);
            color pixel_color(ray_color(r));
            write_color(cout, pixel_color);
        }
    }
    cerr<<endl<<"Done."<<endl;
}
