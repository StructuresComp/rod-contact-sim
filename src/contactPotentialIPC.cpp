#include "contactPotentialIPC.h"

contactPotentialIPC::contactPotentialIPC(std::vector<elasticRod *> m_rod, timeStepper &m_stepper,
                                         collisionDetector &m_col_detector, double m_delta)
                                         : baseContactPotential{m_rod, m_stepper, m_col_detector, m_delta}
                                         {

    p2p_input[6] = delta;
    p2p_input[7] = h2;

    e2p_input[9] = delta;
    e2p_input[10] = h2;

    e2e_input[12] = delta;
    e2e_input[13] = h2;

    bool is_imc = false;
    sym_eqs->generateContactPotentialPiecewiseFunctions(is_imc);
}


void contactPotentialIPC::computeFc() {
    for (int i = 0; i < col_detector->num_collisions; i++) {
        idx1 = col_detector->contact_ids(i, 0);
        idx2 = col_detector->contact_ids(i, 1);
        idx3 = col_detector->contact_ids(i, 2);
        idx4 = col_detector->contact_ids(i, 3);

        constraint_type = static_cast<ConstraintType>(col_detector->contact_ids(i, 4));

        prepContactInput();
        contact_gradient.setZero();

        if (constraint_type == ConstraintType::PointToPoint) {
            sym_eqs->E_p2p_ipc_gradient_func.call(p2p_gradient.data(), p2p_input.data());

            // insert gradient and hessian to contact gradient and contact hessian
            contact_gradient(seq(0, 2)) = p2p_gradient(seq(0, 2));
            contact_gradient(seq(6, 8)) = p2p_gradient(seq(3, 5));
        }
        else if (constraint_type == ConstraintType::PointToEdge) {
            sym_eqs->E_e2p_ipc_gradient_func.call(e2p_gradient.data(), e2p_input.data());

            // insert gradient and hessian to contact gradient and contact hessian
            contact_gradient(seq(0, 2)) = e2p_gradient(seq(0, 2));
            contact_gradient(seq(3, 5)) = e2p_gradient(seq(3, 5));
            contact_gradient(seq(6, 8)) = e2p_gradient(seq(6, 8));
        }
        else if (constraint_type == ConstraintType::EdgeToEdge) {
            sym_eqs->E_e2e_ipc_gradient_func.call(e2e_gradient.data(), e2e_input.data());
            contact_gradient = e2e_gradient;
        }

        contact_gradient *= contact_stiffness;

        // add gradient
        for (int e1 = 0; e1 < 3; e1++) {
            stepper->addForce(4 * (idx1 % nv) + e1, contact_gradient[e1], idx1 / nv);
            stepper->addForce(4 * (idx3 % nv) + e1, contact_gradient[e1 + 3], idx3 / nv);
            stepper->addForce(4 * (idx2 % nv) + e1, contact_gradient[e1 + 6], idx2 / nv);
            stepper->addForce(4 * (idx4 % nv) + e1, contact_gradient[e1 + 9], idx4 / nv);
        }
    }
}


void contactPotentialIPC::computeFcJc() {
    for (int i = 0; i < col_detector->num_collisions; i++) {
        idx1 = col_detector->contact_ids(i, 0);
        idx2 = col_detector->contact_ids(i, 1);
        idx3 = col_detector->contact_ids(i, 2);
        idx4 = col_detector->contact_ids(i, 3);

        constraint_type = static_cast<ConstraintType>(col_detector->contact_ids(i, 4));

        prepContactInput();
        contact_gradient.setZero();
        contact_hessian.setZero();

        if (constraint_type == ConstraintType::PointToPoint) {
            sym_eqs->E_p2p_ipc_gradient_func.call(p2p_gradient.data(), p2p_input.data());
            sym_eqs->E_p2p_ipc_hessian_func.call(p2p_hessian.data(), p2p_input.data());

            // insert gradient and hessian to contact gradient and contact hessian
            contact_gradient(seq(0, 2)) = p2p_gradient(seq(0, 2));
            contact_gradient(seq(6, 8)) = p2p_gradient(seq(3, 5));
            //       edge1 edge3 edge2 edge4
            // edge1
            // edge3
            // edge2
            // edge4
            contact_hessian.block<3, 3>(0, 0) = p2p_hessian.block<3, 3>(0, 0);
            contact_hessian.block<3, 3>(0, 6) = p2p_hessian.block<3, 3>(0, 3);
            contact_hessian.block<3, 3>(6, 0) = p2p_hessian.block<3, 3>(3, 0);
            contact_hessian.block<3, 3>(6, 6) = p2p_hessian.block<3, 3>(3, 3);
        } else if (constraint_type == ConstraintType::PointToEdge) {
            sym_eqs->E_e2p_ipc_gradient_func.call(e2p_gradient.data(), e2p_input.data());
            sym_eqs->E_e2p_ipc_hessian_func.call(e2p_hessian.data(), e2p_input.data());

            // insert gradient and hessian to contact gradient and contact hessian
            contact_gradient(seq(0, 2)) = e2p_gradient(seq(0, 2));
            contact_gradient(seq(3, 5)) = e2p_gradient(seq(3, 5));
            contact_gradient(seq(6, 8)) = e2p_gradient(seq(6, 8));

            //       edge1 edge3 edge2 edge4
            // edge1
            // edge3
            // edge2
            // edge4
            contact_hessian.block<9, 9>(0, 0) = e2p_hessian;
        } else if (constraint_type == ConstraintType::EdgeToEdge) {
            sym_eqs->E_e2e_ipc_gradient_func.call(e2e_gradient.data(), e2e_input.data());
            sym_eqs->E_e2e_ipc_hessian_func.call(e2e_hessian.data(), e2e_input.data());

            contact_gradient = e2e_gradient;
            contact_hessian = e2e_hessian;
        }

        contact_gradient *= contact_stiffness;
        contact_hessian *= contact_stiffness;

        int v1 = idx1 % nv;
        int v2 = idx2 % nv;
        int v3 = idx3 % nv;
        int v4 = idx4 % nv;
        int r1 = idx1 / nv;
        int r2 = idx2 / nv;
        int r3 = idx3 / nv;
        int r4 = idx4 / nv;

        // add gradient
        for (int e1 = 0; e1 < 3; e1++) {
            stepper->addForce(4 * v1 + e1, contact_gradient[e1], r1);
            stepper->addForce(4 * v3 + e1, contact_gradient[e1 + 3], r3);
            stepper->addForce(4 * v2 + e1, contact_gradient[e1 + 6], r2);
            stepper->addForce(4 * v4 + e1, contact_gradient[e1 + 9], r4);
        }

        // add hessian
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                // first row
                stepper->addJacobian(4 * v1 + i, 4 * v1 + j, contact_hessian(j, i), r1, r1);
                stepper->addJacobian(4 * v1 + i, 4 * v3 + j, contact_hessian(3 + j, i), r1, r3);
                stepper->addJacobian(4 * v1 + i, 4 * v2 + j, contact_hessian(6 + j, i), r1, r2);
                stepper->addJacobian(4 * v1 + i, 4 * v4 + j, contact_hessian(9 + j, i), r1, r4);

                // second row
                stepper->addJacobian(4 * v3 + i, 4 * v1 + j, contact_hessian(j, 3 + i), r3, r1);
                stepper->addJacobian(4 * v3 + i, 4 * v3 + j, contact_hessian(3 + j, 3 + i), r3, r3);
                stepper->addJacobian(4 * v3 + i, 4 * v2 + j, contact_hessian(6 + j, 3 + i), r3, r2);
                stepper->addJacobian(4 * v3 + i, 4 * v4 + j, contact_hessian(9 + j, 3 + i), r3, r4);

                // third row
                stepper->addJacobian(4 * v2 + i, 4 * v1 + j, contact_hessian(j, 6 + i), r2, r1);
                stepper->addJacobian(4 * v2 + i, 4 * v3 + j, contact_hessian(3 + j, 6 + i), r2, r3);
                stepper->addJacobian(4 * v2 + i, 4 * v2 + j, contact_hessian(6 + j, 6 + i), r2, r2);
                stepper->addJacobian(4 * v2 + i, 4 * v4 + j, contact_hessian(9 + j, 6 + i), r2, r4);

                // forth row
                stepper->addJacobian(4 * v4 + i, 4 * v1 + j, contact_hessian(j, 9 + i), r4, r1);
                stepper->addJacobian(4 * v4 + i, 4 * v3 + j, contact_hessian(3 + j, 9 + i), r4, r3);
                stepper->addJacobian(4 * v4 + i, 4 * v2 + j, contact_hessian(6 + j, 9 + i), r4, r2);
                stepper->addJacobian(4 * v4 + i, 4 * v4 + j, contact_hessian(9 + j, 9 + i), r4, r4);
            }
        }
    }
}



double contactPotentialIPC::computeUpperBound(const VectorXd &p) {
    double alpha = 1.0;
    if (col_detector->num_collisions == 0)
        return alpha;

    for (int i = 0; i < col_detector->candidate_set.size(); i++) {
        idx1 = col_detector->candidate_set[i][0];
        idx2 = col_detector->candidate_set[i][1];
        idx3 = col_detector->candidate_set[i][0] + 1;
        idx4 = col_detector->candidate_set[i][1] + 1;

        Vector3d v1s = rod_vec[idx1 / nv]->getVertex(idx1 % nv);
        Vector3d v1e = rod_vec[idx3 / nv]->getVertex(idx3 % nv);
        Vector3d v2s = rod_vec[idx2 / nv]->getVertex(idx2 % nv);
        Vector3d v2e = rod_vec[idx4 / nv]->getVertex(idx4 % nv);

        Vector3d p1, p2, p3, p4;
        if (rod_vec[idx1 / nv]->isConstrained[4 * (idx1 % nv)] == 1) {
            p1 = Vector3d(0, 0, 0);
        } else {
            p1 = p.segment((idx1 / nv) * rod_vec[idx1 / nv]->uncons + rod_vec[idx1 / nv]->fullToUnconsMap[4 * (idx1 % nv)], 3);
        }

        if (rod_vec[idx3 / nv]->isConstrained[4 * (idx3 % nv)] == 1) {
            p2 = Vector3d(0, 0, 0);
        } else {
            p2 = p.segment((idx3 / nv) * rod_vec[idx3 / nv]->uncons + rod_vec[idx3 / nv]->fullToUnconsMap[4 * (idx3 % nv)], 3);
        }

        if (rod_vec[idx2 / nv]->isConstrained[4 * (idx2 % nv)] == 1) {
            p3 = Vector3d(0, 0, 0);
        } else {
            p3 = p.segment((idx2 / nv) * rod_vec[idx2 / nv]->uncons + rod_vec[idx2 / nv]->fullToUnconsMap[4 * (idx2 % nv)], 3);
        }

        if (rod_vec[idx4 / nv]->isConstrained[4 * (idx4 % nv)] == 1) {
            p4 = Vector3d(0, 0, 0);
        } else {
            p4 = p.segment((idx4 / nv) * rod_vec[idx4 / nv]->uncons + rod_vec[idx4 / nv]->fullToUnconsMap[4 * (idx4 % nv)], 3);
        }

        double toi;
        bool inContactOrNot = LQF_CCD3(v1s, v1e, v2s, v2e, p1, p2, p3, p4, toi);
        if (inContactOrNot == 1)
            alpha = alpha < toi ? alpha : toi;
    }

    return alpha;
}


bool
contactPotentialIPC::LQF_CCD3(Vector3d &a0s, Vector3d &a1s, Vector3d &b0s, Vector3d &b1s, Vector3d &p1, Vector3d &p2,
                              Vector3d &p3, Vector3d &p4, double &toi, double OrigLen) {
    bool isContact = false;
    int Ntest = 1000;
    double stepsize = 1 / double(Ntest);
    for (int i = 1; i <= Ntest; i++) {
        Vector3d p1_new = i * stepsize * p1;
        Vector3d p2_new = i * stepsize * p2;
        Vector3d p3_new = i * stepsize * p3;
        Vector3d p4_new = i * stepsize * p4;
        Vector3d a0e = a0s + p1_new;
        Vector3d a1e = a1s + p2_new;
        Vector3d b0e = b0s + p3_new;
        Vector3d b1e = b1s + p4_new;
        double dist;
        col_detector->computeMinDistance(a0e, a1e, b0e, b1e, dist);

        // evaluate the results
        if (dist < h2) {
            toi = (i - 1) * stepsize * OrigLen;
            isContact = true;
            if (i == 1) {
                LQF_CCD3(a0s, a1s, b0s, b1s, p1_new, p2_new, p3_new, p4_new, toi, OrigLen / Ntest);
            }
            break;
        }
    }
    return isContact;
}