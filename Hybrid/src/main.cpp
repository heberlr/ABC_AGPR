#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>


#include "Util.hpp"
#include "Cell.hpp"
#include "CellManipulation.hpp"
#include "FileFactory.hpp"
#include "Frame.hpp"
#include "FrameFactory.hpp"
#include "Mesh.hpp"
#include "ConfigHandler.hpp"


#include "Ran.hpp"
#include "NR3.hpp"
#include "Macro.hpp"

int main (int argc, char *argv[]){
	// std::ios::sync_with_stdio(false);
	ConfigHandler config;

	if (!config.created())
		return EXIT_FAILURE;
    
	//One execution
	// Ran ran(config.input.seed);
	//

	//Replica
	Ran ran(time(NULL)+atoi(argv[1]));
    config.output.filenames.number = atoi(argv[1]);
	//

	Frame* frame;
	Mesh* mesh;
	Macro* macro;

	config.parameters.alphaP = atof(argv[2]);
	config.parameters.alphaA = atof(argv[3]);

	//printf("Par: %lf, %lf, %lf, %lf, %lf",config.parameters.alphaP,config.parameters.alphaA,config.agent.oConsumption,config.continuum.oxgD,config.parameters.sigmaH);

	double garbage=0.0;

	if (config.input.fileFormat == BI_D)
	{
		frame = FrameFactory::makeFrame2D(config.input.initialCondition, config.agent.oConsumption); // initial condition
		mesh = new Mesh2D(	frame->domain, config.continuum.oxgD, garbage,
							config.continuum.oBorder, config.continuum.hCoarse,
							config.continuum.hRefined, config.continuum.deltaT	);
		macro = new Macro2D(mesh, frame, &config);
	} else {
		frame = FrameFactory::makeFrame3D(config.input.initialCondition, config.agent.oConsumption); // initial condition
		mesh = new Mesh3D(	frame->domain, config.continuum.oxgD, garbage,
							config.continuum.oBorder, config.continuum.hCoarse,
							config.continuum.hRefined, config.continuum.deltaT	);
		macro = new Macro3D(mesh, frame, &config);
	}

	if (config.output.files)
		FileFactory::makeFrameFile(frame, &config);

	macro->reaction();
	macro->diference();

	//Confluence
/* 	char name[60];
	char path[60];
	sprintf(path,"confluence/confluence%03d/", atoi(argv[2]));
	struct stat info;
	if (stat( path, &info ) != 0) {
		char create[60];
		sprintf(create,"mkdir %s",path);
		system(create);
	}
  sprintf(name,"%s%s%03d-confluence-2D.dat", path,config.output.filenames.agent.c_str(), config.output.filenames.number);
  FILE *arq = fopen(name, "w");
  fprintf(arq, "%05d %e %e\n", frame->time ,macro->LiveConfluence, macro->DeadConfluence); */
    printf("%.1lf %e %e\n", float(frame->time) ,macro->LiveConfluence, macro->DeadConfluence);

	while(frame->time < config.input.timeMax){
		frame->time += config.parameters.delta_tt;

		if (config.input.fileFormat == BI_D)
			CellManipulation2D::updateFrame(frame, &config, mesh, &ran);
		else
			CellManipulation::updateFrame(frame, &config, mesh, &ran);

		if (config.output.prints)
			std::cout << frame->to_string();
		else if(config.output.justLastPrint && frame->time == config.input.timeMax)
			std::cout << frame->to_string();

		if (config.input.fileFormat == BI_D)
			CellManipulation2D::force(frame, &config);
		else
			CellManipulation::force(frame, &config);

		if (config.output.files)
			FileFactory::makeFrameFile(frame, &config);
		else if(config.output.justLastFile && frame->time == config.input.timeMax)
			FileFactory::makeFrameFile(frame, &config);

		macro->reaction();
		macro->diference();
		if (frame->time % 3 == 0){
			//fprintf(arq, "%05d %e %e\n", frame->time, macro->LiveConfluence, macro->DeadConfluence);
            printf("%.1lf %e %e\n", float(frame->time) ,macro->LiveConfluence, macro->DeadConfluence);
        }

	}
	// fclose(arq);
	//Desalocar memória
	delete frame;
	delete mesh;
	delete macro;
	return EXIT_SUCCESS;
}
