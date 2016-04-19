Usage
=====
1 Calculate average model for oceanic plate 
python model_average.py -m newstartingmodel.589.080.JdF2 -o JdF.model -e JdF.model.err  
* m: model from shearwave inversion
* o: output averaged model
* e: output model err for associate model parameter, same structure as in output model

2 Update convert average model to PREM format
python convert_txt2prem.py -i JdF.model -p ocean.4500.iso -o JdF1DQ_0
* i: input model from shearwave inverison or perturbed model for kernel calculation
* p: reference PREM model used for template
* o: updated model in PREM format 
