#include "deploy.h"
#include <iostream>
#include"MCMF.h"
using namespace std;
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    MCMF mm(topo,line_num);
    write_result(mm.NewSolve().c_str(),filename);
}



