#include "baseContactPotential.h"

baseContactPotential::baseContactPotential(vector<elasticRod *> m_rod, timeStepper &m_stepper,
                                           collisionDetector &m_col_detector, double m_delta) {

    rod_vec = m_rod;
    stepper = &m_stepper;
    col_detector = &m_col_detector;
    delta = m_delta;

    h2 = rod_vec[0]->rodRadius * 2;
    nv = rod_vec[0]->nv;
    scale = 1;

    sym_eqs = new symbolicEquations();
}


void baseContactPotential::updateContactStiffness() {
    if (col_detector->candidate_set.size() == 0) return;
    double curr_max_force = 0;
    double curr_force;
    double fx, fy, fz;
    set<int> nodes_to_check;

    // Compute the maximum force that a node experiences.
    for (int i = 0; i < col_detector->candidate_set.size(); i++) {
        nodes_to_check.insert(col_detector->candidate_set[i][0]);
        nodes_to_check.insert(col_detector->candidate_set[i][0] + 1);
        nodes_to_check.insert(col_detector->candidate_set[i][1]);
        nodes_to_check.insert(col_detector->candidate_set[i][1] + 1);
    }

    for (auto i: nodes_to_check) {
        int idx = i % nv;
        int num = i / nv;
        int offset = num * rod_vec[num]->uncons;
        if (rod_vec[num]->getIfConstrained(4 * idx) == 0 &&
            rod_vec[num]->getIfConstrained(4 * idx + 1) == 0 &&
            rod_vec[num]->getIfConstrained(4 * idx + 2) == 0) {
            fx = stepper->Force[offset + rod_vec[num]->fullToUnconsMap[4 * idx]];
            fy = stepper->Force[offset + rod_vec[num]->fullToUnconsMap[4 * idx + 1]];
            fz = stepper->Force[offset + rod_vec[num]->fullToUnconsMap[4 * idx + 2]];
        } else {
            continue;
        }
        curr_force = sqrt(pow(fx, 2) + pow(fy, 2) + pow(fz, 2));
        if (curr_force > curr_max_force) {
            curr_max_force = curr_force;
        }
    }
    contact_stiffness = 1e5 * curr_max_force;
}


void baseContactPotential::prepContactInput() {
    Vector3d x1s = scale * rod_vec[idx1 / nv]->getVertex(idx1 % nv);
    Vector3d x1e = scale * rod_vec[idx3 / nv]->getVertex(idx3 % nv);
    Vector3d x2s = scale * rod_vec[idx2 / nv]->getVertex(idx2 % nv);
    Vector3d x2e = scale * rod_vec[idx4 / nv]->getVertex(idx4 % nv);

    if (constraint_type == ConstraintType::PointToPoint) {
        p2p_input[0] = x1s(0);
        p2p_input[1] = x1s(1);
        p2p_input[2] = x1s(2);

        p2p_input[3] = x2s(0);
        p2p_input[4] = x2s(1);
        p2p_input[5] = x2s(2);
    }

    if (constraint_type == ConstraintType::PointToEdge) {
        e2p_input[0] = x1s(0);
        e2p_input[1] = x1s(1);
        e2p_input[2] = x1s(2);
        e2p_input[3] = x1e(0);
        e2p_input[4] = x1e(1);
        e2p_input[5] = x1e(2);
        e2p_input[6] = x2s(0);
        e2p_input[7] = x2s(1);
        e2p_input[8] = x2s(2);
    }


    if (constraint_type == ConstraintType::EdgeToEdge) {
        e2e_input[0] = x1s(0);
        e2e_input[1] = x1s(1);
        e2e_input[2] = x1s(2);
        e2e_input[3] = x1e(0);
        e2e_input[4] = x1e(1);
        e2e_input[5] = x1e(2);
        e2e_input[6] = x2s(0);
        e2e_input[7] = x2s(1);
        e2e_input[8] = x2s(2);
        e2e_input[9] = x2e(0);
        e2e_input[10] = x2e(1);
        e2e_input[11] = x2e(2);
    }
}

