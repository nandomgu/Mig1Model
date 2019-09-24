function model_g(input)
	function model(dydt, y, parameters, t)

		g=input(t)
	
		k3= parameters[1]
		k2= parameters[2]
		ksnf1= parameters[3]
		ksnf1std1= parameters[4]
		nsnf1= parameters[5]
		nsnf2= parameters[6]
		Snf1tot= parameters[7]
		dmth1= parameters[8]
		nmth1snf3= parameters[9]
		nmth1rgt2= parameters[10]
		dmth1snf3= parameters[11]
		dmth1rgt2= parameters[12]
		smth1= parameters[13]
		kmig1mth1= parameters[14]
		nmig1mth1= parameters[15]
		kmig2mth1= parameters[16]
		nmig2mth1= parameters[17]
		std1tot= parameters[18]
		istd1= parameters[19]
		nstd1= parameters[20]
		estd1max= parameters[21]
		imig1= parameters[22]
		kmig1snf1= parameters[23]
		emig1max= parameters[24]
		dmig2= parameters[25]
		dmig2snf1= parameters[26]
		kmig2snf1= parameters[27]
		smig2= parameters[28]
		kmig2std1= parameters[29]
		nmig2std1= parameters[30]
		kmig2mth1std1= parameters[31]
		nmig2mth1std1= parameters[32]
		dhxt4= parameters[33]
		dhxt4max= parameters[34]
		kdhxt4= parameters[35]
		ndhxt4= parameters[36]
		shxt4= parameters[37]
		khxt4mth1= parameters[38]
		nhxt4mth1= parameters[39]
		khxt4std1= parameters[40]
		nhxt4std1= parameters[41]
		khxt4mth1std1= parameters[42]
		nhxt4mth1std1= parameters[43]
		khxt4mig1= parameters[44]
		khxt4mig2= parameters[45]
		nhxt4mig1= parameters[46]
		nhxt4mig2= parameters[47]
	
		Hxt4= y[1] #init
		Mig1= y[2] #par[48]
		Mig2= y[3] #par[49]
		Mth1= y[4] #par[50]
		Std1= y[5] #par[51]
		y[6]=g #par[53]
	
		Snf1= Snf1tot*(1 + Std1/ksnf1std1)^nsnf1/((1 + Std1/ksnf1std1)^nsnf1 + (g/ksnf1)^nsnf2)
	
		dydt[1]= +(shxt4/(1 + (Mth1/khxt4mth1)^nhxt4mth1 + (Std1/khxt4std1)^nhxt4std1 + (Mth1*Std1/ khxt4mth1std1)^nhxt4mth1std1 + (Mig1/khxt4mig1)^nhxt4mig1 + (Mig2/khxt4mig2)^nhxt4mig2))-dhxt4*Hxt4-(dhxt4max/(1 + (g/kdhxt4)^ndhxt4)*Hxt4)
		dydt[2]= +(imig1*(1 - Mig1))-(emig1max*Snf1*Mig1/(kmig1snf1 + Mig1))
		dydt[3]= +(smig2/(1 + (Mth1/kmig2mth1)^nmig2mth1 + (Std1/kmig2std1)^nmig2std1 + (Mth1*Std1/ kmig2mth1std1)^nmig2mth1std1))-dmig2*Mig2-(dmig2snf1*Snf1*Mig2/(kmig2snf1 + Mig2))
		dydt[4]= +(smth1/(1 + (Mig1/kmig1mth1)^nmig1mth1 + (Mig2/kmig2mth1)^nmig2mth1))-dmth1*Mth1-(dmth1snf3*g^nmth1snf3/(k3^nmth1snf3 + g^nmth1snf3)*Mth1)-(dmth1rgt2*g^nmth1rgt2/(k2^nmth1rgt2 + g^nmth1rgt2)*Mth1)
		dydt[5]= +(istd1*(std1tot - Std1))-(estd1max*g^nstd1/(k2^nstd1 + g^nstd1)*Std1)
		dydt[6]= 0.0
	end
end