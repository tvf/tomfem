#include "tomfem.h"

#include <cassert>
#include <iostream>

namespace {

bool parse_boundary_condition(tomfem::BoundaryCondition *boundary_condition,
                              std::istream& in) {
    size_t v1, v2;
    in >> v1 >> v2;

    boundary_condition->boundary = {--v1, --v2};

    std::string s;
    in >> s;

    if (s == "temperature") {
        float temperature;
        in >> temperature;

        boundary_condition->condition
            = tomfem::DirichletBoundaryCondition{ temperature };

        return true;
    }

    if (s == "flux") {
        float flux;
        in >> flux;

        boundary_condition->condition
            = tomfem::NeumannBoundaryCondition{ flux };

        return true;
    }

    if (s == "convection") {
        std::string keyword;

        in >> keyword;
        if (keyword != "alpha") return false;
        float alpha;
        in >> alpha;

        in >> keyword;
        if (keyword != "ambient") return false;
        float ambient_temperature;
        in >> ambient_temperature;

        boundary_condition->condition
            = tomfem::ConvectionBoundaryCondition{ alpha,
                                                   ambient_temperature };

        return true;
    }

    return false;
}

// Note: parser disaster.
int parse_mesh_and_boundary_conditions(tomfem::Mesh *mesh,
                                       tomfem::BoundaryConditions *bcs,
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
                size_t v1, v2, v3;
                in >> v1 >> v2 >> v3;
                if (in.fail()) return line;
                mesh->elements.push_back({--v1, --v2, --v3});
            } break;

            case '#': {
                std::string s;
                in >> s;
                if (s != "boundary") break;
                in >> s;
                if (s != "condition") break;

                tomfem::BoundaryCondition bc;
                if (!parse_boundary_condition(&bc, in)) return line;
                bcs->push_back(bc);
            } break;

            default: break;
        }

        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        ++line;
    }

    return 0;
}

struct BoundaryConditionFormatter {

    void operator()(const tomfem::DirichletBoundaryCondition& bc) {
        out << "temperature " << bc.temperature;
    };

    void operator()(const tomfem::NeumannBoundaryCondition& bc) {
        out << "flux " << bc.flux;
    };

    void operator()(const tomfem::ConvectionBoundaryCondition& bc) {
        out << "convection alpha " << bc.alpha
            << " ambient " << bc.ambient_temperature;
    };

    std::ostream& out;
};

void write_solution(const tomfem::Mesh& mesh,
                    const tomfem::BoundaryConditions& bcs,
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

    for (size_t i = 0; i < bcs.size(); ++i) {
        out << "# boundary condition "
            << bcs[i].boundary[0] + 1 << ' ' << bcs[i].boundary[1] + 1 << ' ';

        std::visit(BoundaryConditionFormatter{out}, bcs[i].condition);

        out << '\n';
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

    write_solution(mesh, bcs, nodal_temps, std::cout);
}
