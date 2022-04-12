#include "FileFactory.hpp"

void FileFactory::makeFrameFile(Frame* frame, ConfigHandler* config){
    char name[40];
    FILE *arq;

    if (config->input.fileFormat == BI_D) {
        sprintf(name,"%s/%s%d-%05d-2D.dat", config->output.paths.agent.c_str(), config->output.filenames.agent.c_str(), config->output.filenames.number, frame->time);
        //********** Saving the data **********
        arq = fopen(name, "w");
        fprintf(arq, "%lf %lf\n", frame->domain.x, frame->domain.y);
        fprintf(arq, "%ld %d\n", frame->cells.size(), frame->time);
        fprintf(arq, "%d %d\n", frame->outCells, frame->tumorCells);

        for(int i = 0; i < frame->cells.size(); i++){
            fprintf(arq,"%d\n",frame->cells[i].state);
            fprintf(arq,"%lf %lf\n",frame->cells[i].coordinates.x,frame->cells[i].coordinates.y);
            fprintf(arq,"%lf %lf %lf %d %d\n",frame->cells[i].nucleusRadius,frame->cells[i].radius,frame->cells[i].actionRadius,frame->cells[i].lifetime,frame->cells[i].previousState);
            fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].radiusRate,0.0,0.0,frame->cells[i].sigmaO);
            fprintf(arq,"%lf %lf\n", frame->cells[i].speed.x, frame->cells[i].speed.y);
        }
    } else {
        sprintf(name,"%s/%s%d-%05d-3D.dat", config->output.paths.agent.c_str(), config->output.filenames.agent.c_str(), config->output.filenames.number, frame->time);
        //********** Saving the data **********
        arq = fopen(name, "w");
        fprintf(arq, "%lf %lf %lf\n", frame->domain.x, frame->domain.y, frame->domain.z);
        fprintf(arq, "%ld %d\n", frame->cells.size(), frame->time);
        fprintf(arq, "%d %d\n", frame->outCells, frame->tumorCells);

        for(int i = 0; i < frame->cells.size(); i++){
            fprintf(arq,"%d\n",frame->cells[i].state);
            fprintf(arq,"%lf %lf %lf\n",frame->cells[i].coordinates.x,frame->cells[i].coordinates.y, frame->cells[i].coordinates.z);
            fprintf(arq,"%lf %lf %lf %d %d\n",frame->cells[i].nucleusRadius,frame->cells[i].radius,frame->cells[i].actionRadius,frame->cells[i].lifetime,frame->cells[i].previousState);
            fprintf(arq,"%lf %lf %lf %lf %lf\n",frame->cells[i].oConsumption,frame->cells[i].radiusRate,0.0,0.0,frame->cells[i].sigmaO);
            fprintf(arq,"%lf %lf %lf\n", frame->cells[i].speed.x, frame->cells[i].speed.y, frame->cells[i].speed.z);
        }
    }
    fclose(arq);
}

void FileFactory::makeFile(Frame* frame, ConfigHandler* config, Mesh* mesh, OutputMode mode){
    FILE *arq;
    char name[40];

    if (config->input.fileFormat == BI_D) {
        if (mode == NUT) {
            sprintf(name,"%s/%s%d-%05d-2D.dat", config->output.paths.nut.c_str(), config->output.filenames.nut.c_str(), config->output.filenames.number, frame->time);
        }
        arq = fopen(name,"w");
        fprintf(arq,"# %lf %lf\n", frame->domain.x, frame->domain.y);
        fprintf(arq,"# %d %d\n", mesh->unityCoarse.x, mesh->unityCoarse.y);
        fprintf(arq,"# %d %d \n", 0, frame->time);
        if (mode == NUT) {
            for(int i = 0; i < mesh->matrixSize; i++){
                fprintf(arq,"%lf %lf %lf\n",    ( i % (int)mesh->unityCoarse.x ) * mesh->hCoarse,
                                                ( i / (int)mesh->unityCoarse.x ) * mesh->hCoarse,
                                                 mesh->uO[i]);
            }
        }
    } else {
        if (mode == NUT) {
            sprintf(name,"%s/%s%d-%05d-3D.dat", config->output.paths.nut.c_str(), config->output.filenames.nut.c_str(), config->output.filenames.number, frame->time);
        }
        arq = fopen(name,"w");
        fprintf(arq,"# %lf %lf %lf \n", frame->domain.x, frame->domain.y, frame->domain.z);
        fprintf(arq,"# %d %d %d\n", mesh->unityCoarse.x, mesh->unityCoarse.y, mesh->unityCoarse.z);
        fprintf(arq,"# %d %d \n",0, frame->time);
        if (mode == NUT) {
            for(int i = 0; i < mesh->matrixSize; i++){
                fprintf(arq,"%lf %lf %lf %lf\n",    ( ( i % (int)( frame->domain.x * frame->domain.y ) ) % (int)frame->domain.x ) * mesh->hCoarse,
                                                    ( ( i % (int)( frame->domain.x * frame->domain.y ) ) / (int)frame->domain.x ) * mesh->hCoarse,
                                                    ( i / (int)( frame->domain.x * frame->domain.y ) ) * mesh->hCoarse,
                                                    mesh->uO[i]);
            }
        }
    }

    fclose(arq);
}
