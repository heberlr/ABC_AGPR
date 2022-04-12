#include "Macro.hpp"

void Macro3D::diference(){
	double sigma = this->config->continuum.oxgD/pow(this->mesh->hCoarse,2);

	if (config->output.prints) std::cout << "Oxygen Concentration..." << '\n';

	for (int i = 0; i < this->mesh->matrixSize; i++){
		if (((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0 &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse != this->frame->domain.x &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0 &&
			((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse!=this->frame->domain.y &&
			(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse != 0 &&
			(i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse!=this->frame->domain.z)
				set_crs(i, i, &mesh->rowPtr[0], mesh->colInd, mesh->val, this->oUptake[i]+6*sigma);
	}
	// double* a = &this->mesh->uO[0];
	//Calcula U[]
	pmgmres_ilu_cr (this->mesh->matrixSize, this->mesh->nzNum, &mesh->rowPtr[0], mesh->colInd, mesh->val, &this->mesh->uO[0], mesh->B, 1000/*itr_max*/, 5/*mr*/, 1.0e-5/*tol_abs*/, 1.0e-5/*tol_rel*/);

	if(this->config->output.nut)
		FileFactory::makeFile(this->frame, this->config, this->mesh, NUT);

}

void Macro3D::reaction(){
	if (config->output.prints) std::cout << "Reaction/Source term ..." << '\n';

  	std::vector<std::vector<std::vector<double> > > uptakeMesh;

	uptakeMesh.resize(this->mesh->unityRefined.x);
	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		uptakeMesh[i].resize(this->mesh->unityRefined.y);
		for (int j = 0; j < this->mesh->unityRefined.x; ++j){
			uptakeMesh[i][j].resize(this->mesh->unityRefined.z);
		}
	}

	//Calcula termo de captaçãthis->mesh->unityCoarse.z
	int DeadCellCount = 0;
	int LiveCellCount = 0;
	Vector3 R;
	for(int l=0;l<this->frame->cells.size();l++){
		for(int k = ( (this->frame->cells[l].coordinates.z - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;k < (this->frame->cells[l].coordinates.z + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;k++){
			for(int j = ( (this->frame->cells[l].coordinates.y - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;j< (this->frame->cells[l].coordinates.y + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;j++){
				for(int i = ( (this->frame->cells[l].coordinates.x - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;i< (this->frame->cells[l].coordinates.x + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;i++){
					if (k >= this->mesh->unityRefined.z || j >= this->mesh->unityRefined.y || i >= this->mesh->unityRefined.x || k < 0 || i < 0 || j < 0) continue;
					if ((sqrt(pow(this->mesh->refined[i].x-this->frame->cells[l].coordinates.x,2) + pow(this->mesh->refined[j].y-this->frame->cells[l].coordinates.y,2) + pow(this->mesh->refined[k].z-this->frame->cells[l].coordinates.z,2)) <= this->frame->cells[l].radius)){
						R = Vector3 (this->mesh->pos[i].x-this->frame->domain.x/2.0,this->mesh->pos[j].y-this->frame->domain.y/2.0,this->mesh->pos[k].z-this->frame->domain.z/2.0);
						if (uptakeMesh[k][j][i] == 0.0 && (sqrt(pow(R.x,2)+pow(R.y,2)+pow(R.z,2)) <= this->frame->domain.x/2.0)){
							if (this->frame->cells[l].state == 4) DeadCellCount++;
							else LiveCellCount++;
						}
						uptakeMesh[k][j][i] += this->frame->cells[l].oConsumption;

					}
				}
			}
		}
	}

	//Confluence
	this->DeadConfluence = DeadCellCount/(double) this->mesh->MatrixRefinedSize;
	this->LiveConfluence = LiveCellCount/(double) this->mesh->MatrixRefinedSize;

	//Calcula this->mesh->unityCoarse.z termo reaçãthis->mesh->unityCoarse.z na malha do nutriente

	double x, y, z;
	double Uptake_temp,Source_temp;
	for(int i = 0; i < this->mesh->matrixSize; i++){
		Uptake_temp=0.0;
		Source_temp=0.0;
		x = ((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))%this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		y = ((i%(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))/this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		z = (i/(this->mesh->unityCoarse.x*this->mesh->unityCoarse.y))*this->mesh->hCoarse;

		//cout << "Nó = " << i << " PONTO = (" << x << ", " << y << ", " << z << ") ";
		//Nós internos
		if (x != 0.0 && x != this->frame->domain.x && y != 0.0 && y != this->frame->domain.y && z != 0.0 && z != this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2) ; k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2) ; j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(this->mesh->hReason+1,3);
			//cout << "Nós internos!" << endl;
			continue;
		}
		//Nós (---,---,0) ■ Face 1
		if (x != 0.0 && x != this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
							if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
							Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 1!" << endl;
			continue;
		}
		//Nós (0,---,---) ■ Face 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0 && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 2!" << endl;
			continue;
		}
		//Nós (---,0,---) ■ Face 3
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 3!" << endl;
			continue;
		}
		//Nós (this->frame->domain.x,---,---) ■ Face 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 4!" << endl;
			continue;
		}
		//Nós (---,this->frame->domain.y,---) ■ Face 5
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 5!" << endl;
			continue;
		}
		//Nós (---,---,this->frame->domain.z) ■ Face 6
		if (x != 0.0 && x != this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(this->mesh->hReason+1,2)*((this->mesh->hReason/2)+1));
			//cout << "Face 6!" << endl;
			continue;
		}
		//Nós (---,0,0) -- Aresta 1
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 1!" << endl;
			continue;
		}
		//Nós (0,---,0) -- Aresta 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 2!" << endl;
			continue;
		}
		//Nós (---,this->frame->domain.y,0) -- Aresta 3
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 3!" << endl;
			continue;
		}
		//Nós (this->frame->domain.x,---,0) -- Aresta 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 4!" << endl;
			continue;
		}
		//Nós (---,0,this->frame->domain.z) -- Aresta 5
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0 && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 6!" << endl;
			continue;
		}
		//Nós (---,this->frame->domain.y,this->frame->domain.z) -- Aresta 7
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 7!" << endl;
			continue;
		}
		//Nós (this->frame->domain.x,---,this->frame->domain.z) -- Aresta 8
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0 && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 8!" << endl;
			continue;
		}
		//Nós (0,0,---) -- Aresta 9
		if (x == 0.0 && y == 0.0 && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 9!" << endl;
			continue;
		}
		//Nós (0,this->frame->domain.y,---) -- Aresta 10
		if (x == 0.0 && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 10!" << endl;
			continue;
		}
		//Nós (this->frame->domain.x,0,---) -- Aresta 11
		if (x == this->frame->domain.x && y == 0.0 && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 11!" << endl;
			continue;
		}
		//Nós (this->frame->domain.x,this->frame->domain.y,---) -- Aresta 12
		if (x == this->frame->domain.x && y == this->frame->domain.y && z != this->frame->domain.z && z != 0.0){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/(pow(((this->mesh->hReason/2))+1,2)*(this->mesh->hReason+1));
			//cout << "Aresta 12!" << endl;
			continue;
		}
		//Nó (0,0,0) * Ponto 1
		if (x == 0.0 && y == 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = 0;j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 1!" << endl;
			continue;
		}
		//Nó (this->frame->domain.x,0,0) * Ponto 2
		if (x == this->frame->domain.x && y == 0.0 && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 2!" << endl;
			continue;
		}
		//Nó (0,this->frame->domain.y,0) * Ponto 3
		if (x == 0.0 && y == this->frame->domain.y && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 3!" << endl;
			continue;
		}
		//Nó (this->frame->domain.x,this->frame->domain.y,0) * Ponto 4
		if (x == this->frame->domain.x && y == this->frame->domain.y && z == 0.0){
			for(int k = 0; k <= (this->mesh->hReason/2); k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 5!" << endl;
			continue;
		}
		//Nó (this->frame->domain.x,0,this->frame->domain.z) * Ponto 6
		if (x == this->frame->domain.x && y == 0.0 && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 6!" << endl;
			continue;
		}
		//Nó (0,this->frame->domain.y,this->frame->domain.z) * Ponto 7
		if (x == 0.0 && y == this->frame->domain.y && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 7!" << endl;
			continue;
		}
		//Nó (this->frame->domain.x,this->frame->domain.y,this->frame->domain.z) * Ponto 8
		if (x == this->frame->domain.x && y == this->frame->domain.y && z == this->frame->domain.z){
			for(int k = (z/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); k <= (z/this->mesh->hCoarse)*this->mesh->hReason; k++){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[k][l][j] == 0) uptakeMesh[k][l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[k][l][j];
					}
				}
			}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,3);
			//cout << "Ponto 8!" << endl;
			continue;
		}
	}
}


void Macro2D::diference(){
	double sigma = this->mesh->deltaT*this->config->continuum.oxgD/pow(this->mesh->hCoarse,2);
	//Initial Condition
	if (this->frame->time == 0)
		this->mesh->uO = std::vector<double>(this->mesh->matrixSize, this->config->continuum.oBorder);

	if (config->output.prints) std::cout << "Oxygen Concentration..." << '\n';

	for (int i = 0; i < this->mesh->matrixSize; i++){
		if ( ((i%this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0) &&
			 ((i%this->mesh->unityCoarse.x)*this->mesh->hCoarse != this->frame->domain.x) &&
			 ((i/this->mesh->unityCoarse.x)*this->mesh->hCoarse != 0) &&
			 ((i/this->mesh->unityCoarse.x)*this->mesh->hCoarse != this->frame->domain.y) ){
				this->mesh->B[i] = this->mesh->uO[i];
				set_crs(i, i, &mesh->rowPtr[0], mesh->colInd, mesh->val, this->mesh->deltaT*this->oUptake[i]+4*sigma+1.0);
				//if (this->frame->time == 0)  std::cout << "B = " << this->mesh->uO[i] << " border = " << this->config->continuum.oBorder << std::endl;
		}
	}
	// double* a = &this->mesh->uO[0];
	//Calcula U[]
	pmgmres_ilu_cr (this->mesh->matrixSize, this->mesh->nzNum, &mesh->rowPtr[0], mesh->colInd, mesh->val, &this->mesh->uO[0], mesh->B, 1000/*itr_max*/, 5/*mr*/, 1.0e-5/*tol_abs*/, 1.0e-5/*tol_rel*/);

	if(this->config->output.nut)
		FileFactory::makeFile(this->frame, this->config, this->mesh, NUT);
}

void Macro2D::reaction(){
	if (config->output.prints) std::cout << "Reaction/Source term ..." << '\n';

  	std::vector<std::vector<double> > uptakeMesh;

	uptakeMesh.resize(this->mesh->unityRefined.x);
	for (int i = 0; i < this->mesh->unityRefined.x; ++i) {
		uptakeMesh[i].resize(this->mesh->unityRefined.y);
	}

	//Calcula termo de captaçãthis->mesh->unityCoarse.z
	Vector3 R;
	int DeadCellCount = 0;
  int LiveCellCount = 0;
	for(int l=0;l<this->frame->cells.size();l++){
		for(int j = ( (this->frame->cells[l].coordinates.y - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;j< (this->frame->cells[l].coordinates.y + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;j++){
			for(int i = ( (this->frame->cells[l].coordinates.x - this->frame->cells[l].actionRadius)/this->mesh->hRefined) ;i< (this->frame->cells[l].coordinates.x + this->frame->cells[l].actionRadius)/this->mesh->hRefined + 1;i++){
				if (j >= this->mesh->unityRefined.y || i >= this->mesh->unityRefined.x || i < 0 || j < 0) continue;
				if ((sqrt(pow(this->mesh->refined[i].x-this->frame->cells[l].coordinates.x,2) + pow(this->mesh->refined[j].y-this->frame->cells[l].coordinates.y,2) ) <= this->frame->cells[l].radius)){
					R = Vector3 (this->mesh->pos[i].x-this->frame->domain.x/2.0,this->mesh->pos[j].y-this->frame->domain.y/2.0,0.0);
					if (uptakeMesh[j][i] == 0.0 && (sqrt(pow(R.x,2)+pow(R.y,2)+pow(R.z,2)) <= this->frame->domain.x/2.0)){
						if (this->frame->cells[l].state == 4) DeadCellCount++;
						else LiveCellCount++;
					}
					uptakeMesh[j][i] += this->frame->cells[l].oConsumption;
				}
			}
		}
	}

	//Confluence
	this->DeadConfluence = DeadCellCount/(double) this->mesh->MatrixRefinedSize;
	this->LiveConfluence = LiveCellCount/(double) this->mesh->MatrixRefinedSize;

	//Calcula this->mesh->unityCoarse.z termo reaçãthis->mesh->unityCoarse.z na malha do nutriente
	double x, y;
	double Uptake_temp,Source_temp;
	for(int i = 0; i < this->mesh->matrixSize; i++){
		Uptake_temp=0.0;
		Source_temp=0.0;
		x = (i%this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		y = (i/this->mesh->unityCoarse.x)*this->mesh->hCoarse;
		//std::cout << "Nó = " << i << " PONTO = (" << x << ", " << y << ") ";
		//Nós internos
		if (x != 0.0 && x != this->frame->domain.x && y != 0.0 && y != this->frame->domain.y){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2) ; j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/pow(this->mesh->hReason+1,2);
			//std::cout << "Nós internos!" << std::endl;
			continue;
		}
		//Nós (---,0) -- Aresta 1
		if (x != 0.0 && x != this->frame->domain.x && y == 0.0){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/(((this->mesh->hReason/2)+1)*(this->mesh->hReason+1));
			//std::cout << "Aresta 1!" << std::endl;
			continue;
		}
		//Nós (0,---) -- Aresta 2
		if (x == 0.0 && y != this->frame->domain.y && y != 0.0){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/(((this->mesh->hReason/2)+1)*(this->mesh->hReason+1));
			//std::cout << "Aresta 2!" << std::endl;
			continue;
		}
		//Nós (---,this->frame->domain.y) -- Aresta 3
		if (x != 0.0 && x != this->frame->domain.x && y == this->frame->domain.y){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/(((this->mesh->hReason/2)+1)*(this->mesh->hReason+1));
			//std::cout << "Aresta 3!" << std::endl;
			continue;
		}
		//Nós (this->frame->domain.x,---) -- Aresta 4
		if (x == this->frame->domain.x && y != this->frame->domain.y && y != 0.0){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason + (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/(((this->mesh->hReason/2)+1)*(this->mesh->hReason+1));
			//std::cout << "Aresta 4!" << std::endl;
			continue;
		}
		//Nó (0,0) * Ponto 1
		if (x == 0.0 && y == 0.0){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = 0;j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,2);
			//std::cout << "Ponto 1!" << std::endl;
			continue;
		}
		//Nó (this->frame->domain.x,0) * Ponto 2
		if (x == this->frame->domain.x && y == 0.0){
				for(int l = 0; l <= (this->mesh->hReason/2); l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,2);
			//std::cout << "Ponto 2!" << std::endl;
			continue;
		}
		//Nó (0,this->frame->domain.y) * Ponto 3
		if (x == 0.0 && y == this->frame->domain.y){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = 0; j <= (this->mesh->hReason/2); j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,2);
			//std::cout << "Ponto 3!" << std::endl;
			continue;
		}
		//Nó (this->frame->domain.x,this->frame->domain.y) * Ponto 4
		if (x == this->frame->domain.x && y == this->frame->domain.y){
				for(int l = (y/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); l <= (y/this->mesh->hCoarse)*this->mesh->hReason; l++){
					for(int j = (x/this->mesh->hCoarse)*this->mesh->hReason - (this->mesh->hReason/2); j <= (x/this->mesh->hCoarse)*this->mesh->hReason; j++){
						if (uptakeMesh[l][j] == 0) uptakeMesh[l][j] = this->config->continuum.oConsumptionBg;
						Uptake_temp += uptakeMesh[l][j];
					}
				}
			this->oUptake[i] = Uptake_temp/pow(((this->mesh->hReason/2))+1,2);
			//std::cout << "Ponto 5!" << std::endl;
			continue;
		}
	}
}
