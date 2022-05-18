#include "collisionDetector.h"


collisionDetector::collisionDetector(std::vector<elasticRod *> m_rod_vec, timeStepper &m_stepper,
                                     double m_delta, double m_col_limit) {
    rod_vec = m_rod_vec;
    stepper = &m_stepper;
    rod_num = rod_vec.size();
    rod = rod_vec[0];
    nv = rod->nv;

    delta = m_delta;
    col_limit = m_col_limit;
    scale = 1 / rod->rodRadius;

    contact_limit = scale * (2 * rod->rodRadius + delta);
    candidate_limit = scale * (2 * rod->rodRadius + col_limit);
    numerical_limit = scale * (2 * rod->rodRadius - delta);

    num_edge_combos = 0;
    int ignore_adjacent = 3;  // Here, we assume that no edge can collide with the next 3 adjacent edges on either side

    for (int i = 0; i < rod_vec.size() * rod->nv; i++) {
        for (int j = i + 1; j < rod_vec.size() * rod->nv; j++) {
            if (((i + 1) % rod->nv == 0) || ((j + 1) % rod->nv == 0) || j - i < ignore_adjacent + 1) continue;
            num_edge_combos++;
        }
    }

    contact_ids.resize(num_edge_combos, 6);
    edge_ids.resize(num_edge_combos, 2);

    int real_index = 0;
    for (int i = 0; i < rod_vec.size() * rod->nv; i++) {
        for (int j = i + 1; j < rod_vec.size() * rod->nv; j++) {
            if (((i + 1) % rod->nv == 0) || ((j + 1) % rod->nv == 0) || j - i < ignore_adjacent + 1) continue;
            edge_ids(real_index, 0) = i;
            edge_ids(real_index, 1) = j;
            real_index++;
        }
    }
}

void collisionDetector::fixbound(double &x) {
    if (x > 1) {
        x = 1;
    } else if (x < 0) {
        x = 0;
    }
}


void collisionDetector::computeMinDistance(const Vector3d &v1s, const Vector3d &v1e, const Vector3d &v2s,
                                           const Vector3d &v2e, double &dist) {
    Vector3d e1 = v1e - v1s;
    Vector3d e2 = v2e - v2s;
    Vector3d e12 = v2s - v1s;

    double D1 = e1.array().pow(2).sum();
    double D2 = e2.array().pow(2).sum();
    double R = (e1.array() * e2.array()).sum();
    double S1 = (e1.array() * e12.array()).sum();
    double S2 = (e2.array() * e12.array()).sum();

    double den = D1 * D2 - pow(R, 2);

    double t = 0.0;
    if (den != 0) {
        t = (S1 * D2 - S2 * R) / den;
    }
    fixbound(t);

    double u = (t * R - S2) / D2;

    double uf = u;
    fixbound(uf);

    if (uf != u) {
        t = (uf * R + S1) / D1;
    }
    fixbound(t);

    dist = (e1 * t - e2 * uf - e12).norm();
}


void collisionDetector::computeMinDistance(int &idx1, int &idx2, int &idx3, int &idx4, double &dist, ConstraintType &constraint_type) {
    // Here, we use Lumelsky's min-distance algorithm to figure out the constraint type for each contact pair

    Vector3d v1s = rod_vec[idx1 / nv]->getVertex(idx1 % nv) * scale;
    Vector3d v1e = rod_vec[idx1 / nv]->getVertex(idx1 % nv + 1) * scale;

    Vector3d v2s = rod_vec[idx2 / nv]->getVertex(idx2 % nv) * scale;
    Vector3d v2e = rod_vec[idx2 / nv]->getVertex(idx2 % nv + 1) * scale;

    Vector3d e1 = v1e - v1s;
    Vector3d e2 = v2e - v2s;
    Vector3d e12 = v2s - v1s;

    double D1 = e1.array().pow(2).sum();
    double D2 = e2.array().pow(2).sum();
    double R = (e1.array() * e2.array()).sum();
    double S1 = (e1.array() * e12.array()).sum();
    double S2 = (e2.array() * e12.array()).sum();

    double den = D1 * D2 - pow(R, 2);

    double t = 0.0;
    if (den != 0) {
        t = (S1 * D2 - S2 * R) / den;
    }
    fixbound(t);

    double u = (t * R - S2) / D2;

    double uf = u;
    fixbound(uf);

    if (uf != u) {
        t = (uf * R + S1) / D1;
    }
    fixbound(t);

    // arrange the idx of the nodes based on the contact type
    // p2p, idx1 and idx2 are for two contact nodes
    // p2e, idx1 is for p; idx2 and idx3 are for e
    // e2e, idx1 and idx2 are for one e; idx3 and idx4 are for another e
    if ((t == 0 || t == 1) && (uf == 0 || uf == 1))
    {
        if (t == 0) {
            idx3 = idx1 + 1;
        }
        if (t == 1) {
            idx3 = idx1;
            idx1 = idx1 + 1;
        }
        if (uf == 0) {
            idx4 = idx2 + 1;
        }
        if (uf == 1) {
            idx4 = idx2;
            idx2 = idx2 + 1;
        }
        constraint_type = ConstraintType::PointToPoint;
    } else
    {
        if (t == 0 || t == 1 || uf == 0 || uf == 1)
        {
            if (t == 0) {
                int temp = idx1;
                idx1 = idx2;
                idx2 = temp;
                idx3 = idx1 + 1;
                idx4 = idx2 + 1;
            }
            if (t == 1) {
                int temp = idx1 + 1;
                idx1 = idx2;
                idx2 = temp;
                idx3 = idx1 + 1;
                idx4 = idx2 - 1;
            }
            if (uf == 0) {
                idx2 = idx2;
                idx3 = idx1 + 1;
                idx4 = idx2 + 1;
            }
            if (uf == 1) {
                idx2 = idx2 + 1;
                idx3 = idx1 + 1;
                idx4 = idx2 - 1;
            }
            constraint_type = ConstraintType::PointToEdge;
        }
        else {
            idx3 = idx1 + 1;
            idx4 = idx2 + 1;
            constraint_type = ConstraintType::EdgeToEdge;
        }

    }
    dist = (e1 * t - e2 * uf - e12).norm();
}


void collisionDetector::constructCandidateSet() {
    int edge1, edge2;
    double curr_dist;
    min_dist = 1e10;  // something arbitrarily large
    candidate_set.clear();

    for (int i = 0; i < num_edge_combos; i++) {
        edge1 = edge_ids(i, 0);
        edge2 = edge_ids(i, 1);
        computeMinDistance(rod_vec[edge1 / nv]->getVertex(edge1 % nv) * scale,
                           rod_vec[edge1 / nv]->getVertex(edge1 % nv + 1) * scale,
                           rod_vec[edge2 / nv]->getVertex(edge2 % nv) * scale,
                           rod_vec[edge2 / nv]->getVertex(edge2 % nv + 1) * scale, curr_dist);

        if (curr_dist < min_dist) {
            min_dist = curr_dist;
        }
        if (curr_dist < candidate_limit) {
            candidate_set.push_back(Vector2i(edge1, edge2));
        }
    }
    min_dist /= scale;
}


void collisionDetector::detectCollisions() {
    int idx1, idx2, idx3, idx4;
    ConstraintType constraint_type;
    double curr_dist;
    int j = 0;

    for (int i = 0; i < candidate_set.size(); i++) {
        idx1 = candidate_set[i][0];
        idx2 = candidate_set[i][1];
        computeMinDistance(idx1, idx2, idx3, idx4, curr_dist, constraint_type);

        if (curr_dist < contact_limit) {
            contact_ids(j, 0) = idx1;
            contact_ids(j, 1) = idx2;
            contact_ids(j, 2) = idx3;
            contact_ids(j, 3) = idx4;
            contact_ids(j, 4) = constraint_type;
            contact_ids(j, 5) = (curr_dist > numerical_limit) ? ContactPiecewise::Normal : ContactPiecewise::Penetrated;
            j++;
        }

    }
    num_collisions = j;
}
