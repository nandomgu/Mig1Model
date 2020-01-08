using Dierckx
using DifferentialEquations
function wt(params)
   	params
end

include("genotypes.jl")
#Model 6. bonus. removing std1- snf1 interaction and mig -mth1 interaction and the std1-mig2 interaction. 

function makeproblem(modelfile, inputx, inputy, concentration, t, h0, genotype = wt)
    #making a general spline for every purpose
    input = Spline1D(inputx, inputy * concentration; k = 1)
    #SYSTEM OF ODES WITH A TIME VARYING INPUT
    
    function modelg(input)
    
    
        function model(dydt, y, parameters, t)
    
            g = input(t)
    
    		j=1
            k3= parameters[j]; j+=1;#[1]
            k2= parameters[j]; j+=1;#[2]
            ksnf1= parameters[j]; j+=1;#[3]
            #ksnf1std1= parameters[j]; j+=1;#[4]
            #nsnf1= parameters[j]; j+=1;#[5]
            nsnf2= parameters[j]; j+=1;#[6]
            Snf1tot= parameters[j]; j+=1;#[7]
            dmth1= parameters[j]; j+=1;#[8]
            nmth1snf3= parameters[j]; j+=1;#[9]
            nmth1rgt2= parameters[j]; j+=1;#[10]
            dmth1snf3= parameters[j]; j+=1;#[11]
            dmth1rgt2= parameters[j]; j+=1;#[12]
            smth1= parameters[j]; j+=1;#[13]
            #kmth1mig1= parameters[j]; j+=1;#[14]
            #nmth1mig1= parameters[j]; j+=1;#[15]
            #kmth1mig2= parameters[j]; j+=1;#[16]
            #nmth1mig2= parameters[j]; j+=1;#[17]
            std1tot= parameters[j]; j+=1;#[18]
            istd1= parameters[j]; j+=1;#[19]
            nstd1= parameters[j]; j+=1;#[20]
            nstd3= parameters[j]; j+=1;#[21]
            estd1max= parameters[j]; j+=1;#[22]
            imig1= parameters[j]; j+=1;#[23]
            kmig1snf1= parameters[j]; j+=1;#[24]
            emig1max= parameters[j]; j+=1;#[25]
            dmig2= parameters[j]; j+=1;#[26]
            dmig2snf1= parameters[j]; j+=1;#[27]
            kmig2snf1= parameters[j]; j+=1;#[28]
            smig2= parameters[j]; j+=1;#[29]
            #kmig2std1= parameters[j]; j+=1;#[30]
            #nmig2std1= parameters[j]; j+=1;#[31]
            #kmig2mth1std1= parameters[j]; j+=1;#[32]
            #nmig2mth1std1= parameters[j]; j+=1;#[33]
            kmig2mth1= parameters[j]; j+=1;#[34]
            nmig2mth1= parameters[j]; j+=1;#[35]
            dhxt4= parameters[j]; j+=1;#[36]
            dhxt4max= parameters[j]; j+=1;#[37]
            kdhxt4= parameters[j]; j+=1;#[38]
            ndhxt4= parameters[j]; j+=1;#[39]
            shxt4= parameters[j]; j+=1;#[40]
            khxt4mth1= parameters[j]; j+=1;#[41]
            nhxt4mth1= parameters[j]; j+=1;#[42]
            khxt4std1= parameters[j]; j+=1;#[43]
            nhxt4std1= parameters[j]; j+=1;#[44]
            #khxt4mth1std1= parameters[j]; j+=1;#[45]
            #nhxt4mth1std1= parameters[j]; j+=1;#[46]
            khxt4mig1= parameters[j]; j+=1;#[47]
            khxt4mig2= parameters[j]; j+=1;#[48]
            nhxt4mig1= parameters[j]; j+=1;#[49]
            nhxt4mig2= parameters[j]; j+=1;#[50]

        
            Hxt4 = y[1]
            Mig1 = y[2] 
            Mig2 = y[3] 
            Mth1 = y[4] #parameters[51]
            Std1 = y[5] #parameters[52]
            y[6] = g
    		
    		#(1 + Std1/ksnf1std1)^nsnf1
            Snf1= Snf1tot/(1  + (g/ksnf1)^nsnf2)
			#removed the + (Mth1*Std1/ khxt4mth1std1)^nhxt4mth1std1 
            dydt[1]= +(shxt4/(1 + (Mth1/khxt4mth1)^nhxt4mth1 + (Std1/khxt4std1)^nhxt4std1 + (Mig1/khxt4mig1)^nhxt4mig1 + (Mig2/khxt4mig2)^nhxt4mig2))-dhxt4*Hxt4-(dhxt4max/(1 + (g/kdhxt4)^ndhxt4)*Hxt4)
            dydt[2]= +(imig1*(1 - Mig1))-(emig1max*Snf1*Mig1/(kmig1snf1 + Mig1))
            #removed (Mth1*Std1/ kmig2mth1std1)^nmig2mth1std1)
            # removed (Std1/kmig2std1)^nmig2std1)
            dydt[3]= +(smig2/(1 + (Mth1/kmig2mth1)^nmig2mth1 -dmig2*Mig2-(dmig2snf1*Snf1*Mig2/(kmig2snf1 + Mig2))
            #removed repressive terms (Mig1/kmth1mig1)^nmth1mig1 + (Mig2/kmth1mig2)^nmth1mig2)
            dydt[4]= +smth1 -dmth1*Mth1-(dmth1snf3*g^nmth1snf3/(k3^nmth1snf3 + g^nmth1snf3)*Mth1)-(dmth1rgt2*g^nmth1rgt2/(k2^nmth1rgt2 + g^nmth1rgt2)*Mth1)
            dydt[5]= +(istd1*(std1tot - Std1))-(estd1max*g^nstd1/(k2^nstd1 + g^nstd1)*Std1 + estd1max*g^nstd1/(k3^nstd3 + g^nstd3)*Std1)
            dydt[6]= 0.0
        end
    
    end
    #OUTPUTING A FUNCTION THAT GENERATES A PROBLEM AS A FUNCTION OF ONLY PARAMETERS.
    
    prob(params) = ODEProblem(modelg(input), [h0;0.0;0.0;exp.(genotype(params)[[end-1, end]]);0.0], t, [exp(j) for j in genotype(params)[1:end-2]])
    
end
