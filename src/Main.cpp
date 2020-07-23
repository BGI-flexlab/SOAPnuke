#include <iostream>
#include <string>
#include <fstream>

#include "process_argv.h"
#include "global_parameter.h"
#include "peprocess.h"
#include "seprocess.h"
#include "processHts.h"
#include "processStLFR.h"

int main(int argc,char* argv[]){
	//C_global_variable gv;	//global variable, include statistic information
	C_global_parameter gp;	//global parameter generated by commands
	check_module(argc,argv);	//check module
	global_parameter_initial(argc,argv,gp);	//initial global parameter
	//C_global_parameter gp=C_global_parameter(argc,argv);
	check_parameter(argc,argv,gp);	//check global parameter whether correct
	if(gp.module_name=="filterHts"){
	    processHts new_task(gp);
	    if(new_task.pe){
	        new_task.processPE();
	    }else{
	        new_task.processSE();
	    }
//	    new_task.process();
	}else if(gp.module_name=="filterStLFR") {
        peProcess* new_task=new processStLFR(gp);
        new_task->process();
	}else {
        if (!gp.fq2_path.empty()) {    //PE data
            peProcess new_task(gp);
            new_task.process();
        } else {    //SE data
            seProcess new_task(gp);
            new_task.process();
        }
	}
	return 0;
}
