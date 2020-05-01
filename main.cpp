#include "tomfem.h"

#include <cassert>
#include <iostream>

namespace {

// Note: parser disaster.
int parse_mesh_and_boundary_conditions(tomfem::Mesh *mesh,
                                        tomfem::BoundaryConditions * /* bcs */,
                                        std::istream& in) {
    int line = 1;

    while (in) {

        char c;
        in >> c;

        switch (c) {
            case 'v': {
                          float x, y, z;
                          in >> x >> y >> z;
                          if (in.fail()) return line;
                          mesh->nodes.push_back({x, y});
                      } break;

            case 'f': {
                          int v1, v2, v3;
                          in >> v1 >> v2 >> v3;
                          if (in.fail()) return line;
                          mesh->elements.push_back({--v1, --v2, --v3});
                      } break;

            case '#': {
                          // parse boundary conditions
                      } break;

            default: break;
        }

        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        ++line;
    }

    return 0;
}

void write_solution(const tomfem::Mesh& mesh,
                    const std::vector<float>& nodal_temps,
                    std::ostream& out) {

    assert(nodal_temps.size() == mesh.nodes.size());

    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        out << 'v' << ' '
            << mesh.nodes[i].x << ' '
            << mesh.nodes[i].y << ' '
            << nodal_temps[i] << '\n';
    }

    for (size_t i = 0; i < mesh.elements.size(); ++i) {
        out << 'f' << ' '
            << mesh.elements[i][0] + 1 << ' '
            << mesh.elements[i][1] + 1 << ' '
            << mesh.elements[i][2] + 1 << '\n';
    }
}

}  // unnamed namespace

int main(int, char**) {

    tomfem::Mesh mesh;
    tomfem::BoundaryConditions bcs;

    int line = parse_mesh_and_boundary_conditions(&mesh, &bcs, std::cin);

    if (line) {
        std::cout << "error parsing input around line: " << line << '\n';
        return 1;
    }

    std::vector<float> nodal_temps = tomfem::solve_heat_equation(mesh, bcs);

    write_solution(mesh, nodal_temps, std::cout);
}
