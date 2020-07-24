#include "CellManipulation.hpp"

Cell CellManipulation::divide(Cell *cell, double rand1, double rand2){
    double  teta = 2*M_PI*rand1,
    phi = M_PI*rand2,
    radiusNucleus = cell->nucleusRadius/pow(2,(1/3.0));
    Vector3 pos = cell->coordinates, dist(radiusNucleus*sin(phi)*cos(teta), radiusNucleus*sin(phi)*sin(teta), radiusNucleus*cos(phi));

    cell->coordinates.x = pos.x + dist.x;
    cell->coordinates.y = pos.y + dist.y;
    cell->coordinates.z = pos.z + dist.z;
    cell->radius = cell->radius/pow(2,(1/3.0));
    cell->nucleusRadius = radiusNucleus;
    cell->actionRadius = cell->actionRadius/pow(2,(1/3.0));
    cell->state = 5;

    Cell child = *cell;
    child.coordinates.x = pos.x - dist.x;
    child.coordinates.y = pos.y - dist.y;
    child.coordinates.z = pos.z - dist.z;
    return child;
}

void CellManipulation::force(Frame *frame, ConfigHandler *config){
    double iter = 0, domainRadius = frame->domain.y/2;
    if (config->output.prints)
        printf("Forces...");
    int count=0;
    while (iter < 60){
        int n =1; int m = 1; double M = 1;
        for(int i = 0; i < frame->cells.size(); i++){
            for(int j = i+1; j < frame->cells.size(); j++){
                Vector3 r(
                            frame->cells[j].coordinates.x - frame->cells[i].coordinates.x,
                            frame->cells[j].coordinates.y - frame->cells[i].coordinates.y,
                            frame->cells[j].coordinates.z - frame->cells[i].coordinates.z
                        );
                double actionRadius = frame->cells[j].actionRadius + frame->cells[i].actionRadius;

                if (norma(r) > actionRadius) continue;

                Vector3 phi = func_var_phi(r, actionRadius, n),
                        psi = func_var_psi( r,
                                            frame->cells[j].nucleusRadius + frame->cells[i].nucleusRadius,
                                            frame->cells[j].radius + frame->cells[i].radius,
                                            M, m);

                Vector3 f((phi*config->forces.c_cca) + (psi*config->forces.c_ccr));

                frame->cells[i].force += f;
                frame->cells[j].force -= f;
            }

            if ( norma(frame->cells[i].coordinates - domainRadius) < (domainRadius - frame->cells[i].actionRadius) ) continue;
            Vector3 _normal = normal(frame->cells[i].coordinates, domainRadius);
            Vector3 phi = func_var_phi(_normal, frame->cells[i].actionRadius, n),
                    psi = func_var_psi(_normal, frame->cells[i].nucleusRadius, frame->cells[i].radius, M, m);

            Vector3 f((phi*config->forces.c_ct) + (psi*config->forces.c_rct));

            frame->cells[i].force += f;
        }

        //Calcula velocidade
        #pragma omp parallel for
        for(int i = 0; i < frame->cells.size(); i++){
            frame->cells[i].speed.x = -0.5*(frame->cells[i].force.x);
            frame->cells[i].speed.y = -0.5*(frame->cells[i].force.y);
            frame->cells[i].speed.z = -0.5*(frame->cells[i].force.z);
        }

        //Norma da velocidade maior
        double v_max= 10e-10;
        #pragma omp parallel for reduction(max : v_max)
        for(int i = 0; i < frame->cells.size(); i++){
            if (norma(frame->cells[i].speed) >= v_max) v_max = norma(frame->cells[i].speed);
        }
        double delta_tt = (1/v_max);

        #pragma omp parallel for
        for(int i = 0; i < frame->cells.size(); i++){
            //Desloca Células
            frame->cells[i].coordinates.x += delta_tt * frame->cells[i].speed.x;
            frame->cells[i].coordinates.y += delta_tt * frame->cells[i].speed.y;
            frame->cells[i].coordinates.z += delta_tt * frame->cells[i].speed.z;

            //ZeactionRadius somatório de forças
            frame->cells[i].force = Vector3();

        }

        for(int i = 0; i < frame->cells.size(); i++){
            double dist = sqrt( pow(frame->cells[i].coordinates.x-domainRadius,2)+pow(frame->cells[i].coordinates.y-domainRadius,2)+pow(frame->cells[i].coordinates.z-domainRadius,2) );
            //Se passar do domínio exclui
            if (dist >= domainRadius){
                if (frame->cells[i].state != 6) frame->tumorCells -= 1;
                frame->cells.erase (frame->cells.begin() + i);
                i--;
                frame->outCells += 1;
            }
        }

        iter+= delta_tt;
        count++;
    }
    if (config->output.prints)
	{
        printf(" (%d iterations)\n",count);
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Norma euclidiana
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline double CellManipulation::norma(Vector3 pos){
    return sqrt(pow(pos.x,2) + pow(pos.y,2) + pow(pos.z,2));
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Calcula vetor normal
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation::normal(Vector3 coordinates, double domainRadius){
    //Direção
    Vector3 u(coordinates.x - domainRadius, coordinates.y - domainRadius, coordinates.z - domainRadius);
    // double u1 = , u2 =  u3 =;
    double D = norma(u);
    //Vesor
    u /= D; // u1 = u1/D; u2 = u2/D; u3 = u3/D;
    //Ponto na casca esférica
    u *= domainRadius; // double X = u1*domainRadius, Y = u2*domainRadius, Z = u3*domainRadius;

    return coordinates - u; //Vector3(x_cel - u.x, y_cel - u.y, z_cel - u.z);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Função phi
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation::func_var_phi(Vector3 normal, double actionRadius, int n){
    double r = norma(normal);
    if ((r > 0) && (r <= actionRadius)){
        double var = pow((1 - r/actionRadius), n+1)/r;
        // return Vector3(var*normal.x, var*normal.y, var*normal.z);
        return normal * var;
    }

    return Vector3();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Função psi
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vector3 CellManipulation::func_var_psi(Vector3 normal, double nucleusRadius, double radius, double M, int m){
    double r = norma(normal);
    if ((r > 0) && (r < nucleusRadius)){
        double c = pow((1 - nucleusRadius/radius),m+1) - M;
        double var = -(c*r/nucleusRadius + M)/r;
        return Vector3(var*normal.x, var*normal.y, var*normal.z);
    }
    if ((r >= nucleusRadius) && (r <= radius)){
        double var = -pow((1- (r/radius)),m+1)/r;
        return Vector3(var*normal.x, var*normal.y, var*normal.z);
    }
    // if ( (r > radius) || (r == 0) )
    return Vector3();
}

void CellManipulation::calculateCellSigmas(Cell *cell, Mesh *mesh){
    Vector3 pos((int)(cell->coordinates.x/mesh->hCoarse),
                (int)(cell->coordinates.y/mesh->hCoarse),
                (int)(cell->coordinates.z/mesh->hCoarse));
    Vector3 c   (((cell->coordinates.x - ((int)cell->coordinates.x))/mesh->hCoarse),
                ((cell->coordinates.y - ((int)cell->coordinates.y))/mesh->hCoarse),
                ((cell->coordinates.z - ((int)cell->coordinates.z))/mesh->hCoarse));

    double oSigma1 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                (int)pos.y*mesh->unityCoarse.y + (int)pos.x];
    double oSigma2 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                (int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
    double oSigma3 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                ((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
    double oSigma4 = mesh->uO[  (int)pos.z*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                ((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];
    double oSigma5 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                (int)pos.y*mesh->unityCoarse.y + (int)pos.x];
    double oSigma6 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                (int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
    double oSigma7 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                ((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
    double oSigma8 = mesh->uO[  ((int)pos.z+1)*(mesh->unityCoarse.y*mesh->unityCoarse.x) +
                                ((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];

    cell->sigmaO = (1-c.x)*(1-c.y)*(1-c.z)*oSigma1 + c.x*(1-c.y)*(1-c.z)*oSigma2 + (1-c.x)*c.y*(1-c.z)*oSigma3 + c.x*c.y*(1-c.z)*oSigma4 + (1-c.x)*(1-c.y)*c.z*oSigma5 + c.x*(1-c.y)*c.z*oSigma6 + (1-c.x)*c.y*c.z*oSigma7 + c.x*c.y*c.z*oSigma8;
}


void CellManipulation::updateFrame(Frame *frame, ConfigHandler *config, Mesh *mesh, Ran *ran)
{
    std::vector<Cell> cellsInserted;
    std::vector<double> randomNumbers = std::vector<double>(frame->cells.size()+1, 0.0);
    // double randomNumbers[frame->cells.size()+1];
    for(int i = 0; i < frame->cells.size()+1; i++)
        randomNumbers[i] = ran->doub();
    #pragma omp parallel for schedule(static) ordered
    for(int i = 0; i < frame->cells.size(); i++){
        calculateCellSigmas(&frame->cells[i], mesh);
        if (frame->cells[i].state == 1){
            //Transiçãmesh.unityCoarse.z Quiescente para Proliferativa
            double alpha_p = config->parameters.alphaP*(frame->cells[i].sigmaO-config->parameters.sigmaH)/(1 - config->parameters.sigmaH);
            alpha_p = max(alpha_p,0.0);
	          double alpha_a = config->parameters.alphaA + config->parameters.gammaA/(1.0+exp(100*(frame->cells[i].sigmaO-config->parameters.sigmaH)));
            if (1-exp(-config->parameters.delta_tt*alpha_p) >= randomNumbers[i]){
                frame->cells[i].lifetime = frame->time;
                frame->cells[i].previousState = frame->cells[i].state;
                frame->cells[i].state = 2;
            }
            //Transiçãmesh.unityCoarse.z Quiescente para Apoptótica
            else if ((1-exp(-config->parameters.delta_tt*alpha_a) >= randomNumbers[i])){
                    frame->cells[i].lifetime = frame->time;
                    frame->cells[i].previousState = frame->cells[i].state;
                    frame->cells[i].oConsumption = 0.0;
                    frame->cells[i].state = 4;
            }
            continue;
        }

        double tau = frame->time - frame->cells[i].lifetime;

        //Ocorre a divisão celular
        if (frame->cells[i].state == 2 && tau == 9){
            //frame->cells.push_back(this->divide(&frame->cells[i], config->input.seed));
            #pragma omp ordered
                cellsInserted.push_back(divide(&frame->cells[i], randomNumbers[i],randomNumbers[i+1]));
            frame->tumorCells++;
            continue;
        }

        //Ocorre a G1
        if ( frame->cells[i].state == 5 && tau > (config->parameters.tauP-config->parameters.tauG1) ){
            if(tau < config->parameters.tauP){
                frame->cells[i].nucleusRadius = config->agent.nucleusRadius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
                frame->cells[i].radius        = config->agent.radius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
                frame->cells[i].actionRadius  = frame->cells[i].radius*1.214;
            }
            else{
                frame->cells[i].previousState = 2;
                frame->cells[i].nucleusRadius = config->agent.nucleusRadius;
                frame->cells[i].radius        = config->agent.radius;
                frame->cells[i].actionRadius  = config->agent.actionRadius;
                frame->cells[i].state         = 1;
                frame->cells[i].lifetime      = frame->time;
            }
        }
    }
    //Insere Células
    for(int i = 0; i < cellsInserted.size(); i++){
        frame->cells.push_back(cellsInserted[i]);
    }
}


Vector3 CellManipulation2D::normal(Vector3 coordinates, double domainRadius){
    //Direção
    Vector3 u(coordinates.x - domainRadius, coordinates.y - domainRadius, 0.0);
    // double u1 = , u2 =  u3 =;
    double D = norma(u);
    //Vesor
    u /= D; // u1 = u1/D; u2 = u2/D; u3 = u3/D;
    //Ponto na casca esférica
    u *= domainRadius; // double X = u1*domainRadius, Y = u2*domainRadius, Z = u3*domainRadius;

    return coordinates - u; //Vector3(x_cel - u.x, y_cel - u.y, z_cel - u.z);
}


void CellManipulation2D::force(Frame *frame, ConfigHandler *config){
    double iter = 0, domainRadius = frame->domain.y/2;
    if (config->output.prints)
        printf("Forces...");
    int count=0;
    while (iter < 60){
        int n =1; int m = 1; double M = 1;
        //#pragma omp parallel for
        for(int i = 0; i < frame->cells.size(); i++){
            for(int j = i+1; j < frame->cells.size(); j++){
                double F_cca_i[2]={0,0},F_ccr_i[2]={0,0};
                Vector3 r(
                            frame->cells[j].coordinates.x - frame->cells[i].coordinates.x,
                            frame->cells[j].coordinates.y - frame->cells[i].coordinates.y,
                            0.0
                        );
                double actionRadius = frame->cells[j].actionRadius + frame->cells[i].actionRadius;
                double nucleusRadius = frame->cells[j].nucleusRadius + frame->cells[i].nucleusRadius;
                double radius   = frame->cells[j].radius + frame->cells[i].radius;

                if (norma(r) > actionRadius) continue;

                // double phi_x,phi_y,phi_z,psi_x,psi_y,psi_z;
                Vector3 phi = func_var_phi(r, actionRadius, n),
                        psi = func_var_psi(r, nucleusRadius, radius, M, m);

                // double c_cca;
                // //if (frame->cells[i].state != 6 && frame->cells[j].state != 6)
                // c_cca = 0.488836;
                // //else
                // //c_cca = 0.588836;
                //double c_ccr = 10;

                //Força adesão célula-célula (i)
                F_cca_i[0] = config->forces.c_cca*phi.x;
                F_cca_i[1] = config->forces.c_cca*phi.y;
                //Força repulsão célula-célula (i)
                F_ccr_i[0] = config->forces.c_ccr*psi.x;
                F_ccr_i[1] = config->forces.c_ccr*psi.y;
                //#pragma omp critical
                //{
                //Somatório de forças
                frame->cells[i].force.x += (F_cca_i[0] + F_ccr_i[0]);
                frame->cells[i].force.y += (F_cca_i[1] + F_ccr_i[1]);

                frame->cells[j].force.x -= (F_cca_i[0] + F_ccr_i[0]);
                frame->cells[j].force.y -= (F_cca_i[1] + F_ccr_i[1]);
                //}
            }
            //Força de rigidez do tecido na célula (i)
            double F_ct[2]={0,0},F_rct[2]={0,0};
            // double K = 0.1;
            // double c_ct = 10.0 * K;
            // double c_rct = 4.88836 * K;
            // double phi_x2,phi_y2,phi_z2,psi_x2,psi_y2,psi_z2;

            Vector3 Dif(frame->cells[i].coordinates.x - domainRadius,frame->cells[i].coordinates.y - domainRadius,0.0);

            if ( norma(Dif) < (domainRadius - frame->cells[i].actionRadius) ) continue;
            Vector3 _normal = normal(frame->cells[i].coordinates, domainRadius);
            Vector3 phi2 = func_var_phi(_normal, frame->cells[i].actionRadius, n),
                    psi2 = func_var_psi(_normal, frame->cells[i].nucleusRadius, frame->cells[i].radius, M, m);

            //Adesão
            F_ct[0] = config->forces.c_ct*phi2.x;
            F_ct[1] = config->forces.c_ct*phi2.y;
            //Repulsão
            F_rct[0] = config->forces.c_rct*psi2.x;
            F_rct[1] = config->forces.c_rct*psi2.y;
            //Somatório de forças
            frame->cells[i].force.x += (F_ct[0] + F_rct[0]);
            frame->cells[i].force.y += (F_ct[1] + F_rct[1]);
        }

        //Calcula velocidade
        #pragma omp parallel for
        for(int i = 0; i < frame->cells.size(); i++){
            frame->cells[i].speed.x = -0.5*(frame->cells[i].force.x);
            frame->cells[i].speed.y = -0.5*(frame->cells[i].force.y);
        }

        //Norma da velocidade maior
        double v_max= 10e-10;
        #pragma omp parallel for reduction(max : v_max)
        for(int i = 0; i < frame->cells.size(); i++){
            if (norma(frame->cells[i].speed) >= v_max) v_max = norma(frame->cells[i].speed);
        }
        double delta_tt = (1/v_max);

        #pragma omp parallel for
        for(int i = 0; i < frame->cells.size(); i++){
            //Desloca Células
            frame->cells[i].coordinates.x += delta_tt * frame->cells[i].speed.x;
            frame->cells[i].coordinates.y += delta_tt * frame->cells[i].speed.y;

            //ZeactionRadius somatório de forças
            frame->cells[i].force = Vector3();

        }

        for(int i = 0; i < frame->cells.size(); i++){
            double dist = sqrt( pow(frame->cells[i].coordinates.x-domainRadius,2)+pow(frame->cells[i].coordinates.y-domainRadius,2));
            //Se passar do domínio exclui
            if (dist >= domainRadius){
                if (frame->cells[i].state != 6) frame->tumorCells -= 1;
                frame->cells.erase (frame->cells.begin() + i);
                i--;
                frame->outCells += 1;
            }
        }

        iter+= delta_tt;
        count++;
    }
    if (config->output.prints)
	{
        printf(" (%d iterations)\n",count);
    }
}

Cell CellManipulation2D::divide(Cell *cell, double rand1){
    double  teta = 2*M_PI*rand1,
    radiusNucleus = cell->nucleusRadius/sqrt(2);
    Vector3 pos = cell->coordinates, dist(radiusNucleus*cos(teta), radiusNucleus*sin(teta), 0.0);

    cell->coordinates.x = pos.x + dist.x;
    cell->coordinates.y = pos.y + dist.y;
    cell->radius = cell->radius/sqrt(2);
    cell->nucleusRadius = radiusNucleus;
    cell->actionRadius = cell->actionRadius/sqrt(2);
    cell->state = 5;

    Cell child = *cell;
    child.coordinates.x = pos.x - dist.x;
    child.coordinates.y = pos.y - dist.y;
    return child;
}

void CellManipulation2D::calculateCellSigmas(Cell *cell, Mesh *mesh){
    Vector3 pos((int)(cell->coordinates.x/mesh->hCoarse),
                (int)(cell->coordinates.y/mesh->hCoarse),
                0.0);
    Vector3 c   (((cell->coordinates.x - ((int)cell->coordinates.x))/mesh->hCoarse),
                ((cell->coordinates.y - ((int)cell->coordinates.y))/mesh->hCoarse),
                0.0);

    double oSigma1 = mesh->uO[  (int)pos.y*mesh->unityCoarse.y + (int)pos.x];
    double oSigma2 = mesh->uO[  (int)pos.y*mesh->unityCoarse.y + ((int)pos.x+1)];
    double oSigma3 = mesh->uO[  ((int)pos.y+1)*mesh->unityCoarse.y + (int)pos.x];
    double oSigma4 = mesh->uO[  ((int)pos.y+1)*mesh->unityCoarse.y + ((int)pos.x+1)];

    cell->sigmaO = (1-c.x)*(1-c.y)*oSigma1 + c.x*(1-c.y)*oSigma2 + (1-c.x)*c.y*oSigma3 + c.x*c.y*oSigma4;
}

void CellManipulation2D::updateFrame(Frame *frame, ConfigHandler *config, Mesh *mesh, Ran *ran)
{
    std::vector<Cell> cellsInserted;
    std::vector<double> randomNumbers = std::vector<double>(frame->cells.size()+1, 0.0);
    // double randomNumbers[frame->cells.size()+1];
    for(int i = 0; i < frame->cells.size()+1; i++)
        randomNumbers[i] = ran->doub();
    #pragma omp parallel for schedule(static) ordered
    for(int i = 0; i < frame->cells.size(); i++){
        calculateCellSigmas(&frame->cells[i], mesh);

        if (frame->cells[i].state == 1){
            //Transiçãmesh.unityCoarse.z Quiescente para Proliferativa
            double alpha_p = config->parameters.alphaP*(frame->cells[i].sigmaO-config->parameters.sigmaH)/(1 - config->parameters.sigmaH);
            alpha_p = max(alpha_p,0.0);
	          double alpha_a = config->parameters.alphaA*1.0/(1.0+exp(100*(frame->cells[i].sigmaO-config->parameters.sigmaH)));
            if (1-exp(-config->parameters.delta_tt*alpha_p) >= randomNumbers[i]){
                frame->cells[i].lifetime = frame->time;
                frame->cells[i].previousState = frame->cells[i].state;
                frame->cells[i].state = 2;
                continue;
            }
            //Transiçãmesh.unityCoarse.z Quiescente para Apoptótica
            else if ((1-exp(-config->parameters.delta_tt*alpha_a) >= randomNumbers[i])){
                    frame->cells[i].lifetime = frame->time;
                    frame->cells[i].previousState = frame->cells[i].state;
                    frame->cells[i].oConsumption = 0.0;
                    frame->cells[i].state = 4;
                    continue;
            }
        }

        double tau = frame->time - frame->cells[i].lifetime;

        //Ocorre a divisão celular
        if (frame->cells[i].state == 2 && tau == 9){
            //frame->cells.push_back(divide(&frame->cells[i], config->input.seed));
            #pragma omp ordered
                cellsInserted.push_back(divide(&frame->cells[i], randomNumbers[i]));
            frame->tumorCells++;
            continue;
        }

        //Ocorre a G1
        if ( frame->cells[i].state == 5 && tau > (config->parameters.tauP-config->parameters.tauG1) ){

            if(tau < config->parameters.tauP){
                frame->cells[i].nucleusRadius = config->agent.nucleusRadius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
                frame->cells[i].radius = config->agent.radius*pow(0.5*(1+ (config->parameters.tauG1+(tau-config->parameters.tauP))/config->parameters.tauG1),(1/3.0));
                frame->cells[i].actionRadius = frame->cells[i].radius*1.214;
            }
            else{
                frame->cells[i].previousState = 2;
                frame->cells[i].nucleusRadius = config->agent.nucleusRadius;
                frame->cells[i].radius = config->agent.radius;
                frame->cells[i].actionRadius = config->agent.actionRadius;
                frame->cells[i].state = 1;
                frame->cells[i].lifetime = frame->time;
            }
        }
    }

    //Insere Células
    for(int i = 0; i < cellsInserted.size(); i++){
        frame->cells.push_back(cellsInserted[i]);
    }
}
