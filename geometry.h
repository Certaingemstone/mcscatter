// geometry.h : objects that interact with photons
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Dense>
#include "material.h"

using Eigen::Vector3d; 

class Object
{
    public:
        Object();
        Object(Material mat);
        Material material;
        bool is_inside(Vector3d* point);
};

class Cylinder : public Object
{
    private:
        Vector3d center;
        Vector3d axis; // normalized vector 
        Eigen::Matrix3d rotation_inv; // inverse of matrix that rotated axis from (0,0,1)
        double length;
        double radius2;
    public:
        Cylinder(Vector3d center, Eigen::Matrix3d rotation, double length, double radius, Material mat);
        bool is_inside(Vector3d* point); // if point is inside the cylinder
        // std::optional<Vector3d> ray_intersect(Vector3d &point, Vector3d &direction);
};

inline Object::Object() {
    Material placeholder_mat = Material(10000,1);
    material = placeholder_mat;
}

inline Object::Object(Material mat) {
    material = mat;
}

inline bool Object::is_inside(Vector3d* point) {
    return false;
}

inline Cylinder::Cylinder(Vector3d cent, Eigen::Matrix3d rot, double len, double radius, Material mat) {
    // Rotate axis from z axis using rotation matrix. Then translate to center.
    Vector3d u {{0,0,1}};
    axis = (rot * u).normalized();
    rotation_inv = rot.inverse();
    center = cent;
    length = len;
    radius2 = std::pow(radius, 2);
    material = mat;
}

inline bool Cylinder::is_inside(Vector3d* point) {
    // Check if point is inside cylinder
    Vector3d p = *point - center;
    double proj = axis.dot(p); // projection onto axis
    double distance2 = (p - proj * axis).squaredNorm(); // distance square to axis
    return distance2 < radius2 && std::abs(proj) < (length / 2.0);
}

/*
inline std::optional<Vector3d> Cylinder::ray_intersect(Vector3d &point, Vector3d &direction) {
    // Check if ray with direction from point intersects Cylinder.
    // If it does, return a point outside of Cylinder along direction BEFORE the collision
    Vector3d x2 = center + axis;

    return {};
}
*/

#endif
