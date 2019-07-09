function [lsqdiff, y]=multiInputSimulation( modelname, params,outvar, initialConditions, inputs, outputs)
%% inputs is a cell array where each element is an input signal. the output will be the costs (x) and the simulation outputs (y) of the same model across all inputs. all inputs will be oriented to be in the same direct
modelFeatures= extractModelFeatures(modelname); 

model=eval(strjoin({'@' modelname}, ''));
if nargin<4
    initialConditions= zeros(1, modelFeatures.numInit);
end

if nargin<5
complexRamp= [0;0.0223032063571172;0.0405047314534221;0.0577574155317936;0.0708017171371265;0.074879900620598;0.0758356059661576;0.0771340928584555;0.0770603734213315;0.0766017074310701;0.0825233317190962;0.0909744554418606;0.093993550248502;0.102173022231747;0.116563929130705;0.143203882900438;0.176368881684685;0.223713693914748;0.281052257734873;0.346466307100618;0.407476155958698;0.471069387757387;0.533242734001389;0.593344344903071;0.646896154803748;0.701085931884189;0.738849711878534;0.778092173550209;0.811575003034639;0.840202031942638;0.856290962143873;0.882968428322724;0.897828449866704;0.904189374798246;0.904315649063733;0.906949983540048;0.904302812442446;0.89571826831523;0.884935318073253;0.875509014907067;0.86321464541776;0.846263132839432;0.83026705376725;0.813688756593054;0.79768698434489;0.778942343617509;0.756122945858676;0.728587294009643;0.695907127396961;0.653422645982901;0.598590983814381;0.539549562479458;0.479394232916031;0.414454994302869;0.348601794470689;0.285308136848455;0.224213405996943;0.167090244832356;0.122325009547493;0.0901200542782974;0.0726337968136678;0.0647512577549924;0.0600383055243802;0.0550657572829314;0.0492174819494744;0.0441327496532019;0.0404784360648094;0.0386288373598618;0.0392827393567328;0.0416060245174376;0.04385917090531;0.0492989753395697;0.0573014894896179;0.0739614838095577;0.0994063835142636;0.131546699085503;0.168006119875606;0.20781938996431;0.248857320083113;0.290178889871922;0.335607409038797;0.388906458123698;0.449891255385504;0.511048599375835;0.575290050983856;0.639074832665356;0.701137058482281;0.764278106678239;0.82464351652881;0.87528951712772;0.91648686201527;0.945361227023334;0.961844754355451;0.972233559512459;0.981067057337478;0.98776798745466;0.989072093065479;0.983223292248401;0.975825441497439;0.959369224507547;0.941174416237702;0.927473435570653;0.914943640537118;0.902157327430209;0.893530447953199;0.882762771839404;0.869213756591971;0.858531802932826;0.844644637638763;0.820074164859752;0.790494804327956;0.752437898339734;0.701541856613617;0.635115335626438;0.57337245957311;0.506332572458254;0.437997050295204;0.374797729083169;0.318603618410695;0.263920417707352;0.219355412800032;0.181480990485946;0.147997849252685;0.123091198645273;0.104051523676733;0.0870487941621835;0.0738718773828963;0.0639971956947056;0.0562793281769372;0.0495754666431618;0.0448937810251198;0.0412568388348441;0.0384279360909205;0.0358014676960468;0.0335895169856507;0.0316324820712432;0.0297321354517882;0.0275405856334663;0.0256602563979294;0.0234919909693372;0.0217603007531544;0.0201354702808326;0.0186647892246733;0.017383650443999;0.016511174268384;0.0152562395520462;0.0142146133505696;0.013299036415849;0.0122787686065637;0.0113007857746887;0.0104571487391209;0.00966251216427188;0.00902375582022215;0.00839836325494823;0.00778158911118866;0.00718619652525696]';
stepInput=[0.0320436786099739 0.0119079187717726 0.00543074381957539 0.00355914726328773 0.00243231673854792 8.85984739511592e-05 0.000293693454831295 0.00303965627042515 0.000662171490507482 0 0.00318229303344181 0.00430922458115919 0.00490633765625974 0.00633457859779492 0.00588177391544824 0.00394862204545212 0.00610834754625382 0.00565993007093941 0.00502727642761124 0.00436720344299299 0.00422003734455592 0.00441807776492646 0.00588810724080566 0.00801153809442583 0.00772973929265948 0.00808329044761265 0.00827132589252642 0.00753752184908845 0.00873005387924769 0.0114593189314807 0.0101489580336488 0.0123283800301401 0.0106623645051206 0.0137058831498231 0.011255850431423 0.0134153230652561 0.0138602857242326 0.0136425232646226 0.0170526740006103 0.0216014834623877 0.0247146295580629 0.0295334464471976 0.0324770568373891 0.0333951454894278 0.0407035702833887 0.0493652716330173 0.0756100162020046 0.174282919711837 0.412858248221435 0.620079189695464 0.76536681598833 0.844756910840671 0.899550781577864 0.92183034670076 0.93941200204495 0.947580232752947 0.952776630667009 0.954653476171245 0.95937034634789 0.963255543177662 0.978557975233553 0.982833867739713 0.993728691920712 0.977697456427566 0.99535153159375 0.993445328496454 0.989425510665855 0.976220168715775 0.98430287855909 0.976442359438124 0.998636958615089 0.997654003421186 0.99120683834598 0.982798903798506 0.999525489915838 0.988061484190317 0.992595419863103 1 0.987816934147853 0.994976044721816 0.97832836253113 0.98645792495114 0.992095094510505 0.989540598414687 0.981009549239477 0.987771196315921 0.984781263909241 0.97635487016228 0.963906702356935 0.979786956280874 0.971027652301316 0.973478642717344 0.970431653220771 0.974928562351077 0.981589754436312 0.952639261179128 0.96491987343329 0.966719804395045 0.958634855039702 0.961879050602964 0.968287106518936 0.977077026330968 0.966912210812339 0.951437877261639 0.959139041242801 0.95390781326874 0.963850402649935 0.953104665936031 0.948358992901325 0.956389087034516 0.942594011016547 0.94015744009697 0.956519106457244 0.949571981519532 0.932385079641736 0.955264413666682 0.940089303696556 0.955865349478343 0.938686875464951 0.907327030375258 0.929865469206359 0.944492661076056 0.944355023959091 0.937297749569702 0.946390265535185 0.947437520497974 0.933967049099314 0.92272832656275 0.934901849924432 0.945930481747221 0.925574665151147 0.938714462426072 0.927394148751973 0.917605848724732 0.874461452151652 0.66689489698016 0.416046168534046 0.305238391034 0.254110485693981 0.223483905758883 0.204524682248122 0.188714397999524 0.175669476751635 0.166358153307421 0.159325811503011 0.151075126653098 0.145477903458097 0.139483361339639 0.135982018391916 0.133028493088032 0.126091576940796 0.12416019541886 0.117360519743492 0.115224301153959 0.111133794308671 0.109069743270181 0.107675523811607 0.105690722619235 0.10503554520979 0.0997346601851089 0.0980088195656689 0.0970275513341154 0.0930982766142305 0.0926803663399782 0.0901105570800065 0.0910653079264758 0.0841816516633323 0.0844044389049042 0.0842280952426661 0.0825800587759985 0.0808721529508442 0.0790584118410701 0.0735263614676546 0.0764828070379757 0.0776747645273045 0.0749691111446409 0.0763423666515299 0.0737586241247536 0.0716302577849381 0.0710390016431943];
rampUp=[0 0.00244034094350791 0.0243068310104219 0.0554379390981066 0.0655212567368327 0.0619651318340717 0.060433528996081 0.0554453669215559 0.0514556312723031 0.050885039722151 0.0496655233348287 0.0471682478172535 0.0479595769885814 0.0491847854373926 0.0486142891329156 0.0510104810334834 0.048958856927114 0.0471681145832048 0.0473226104461677 0.0486288673371942 0.0521517495486088 0.0477111836654159 0.0497687787587974 0.0595873830033663 0.0732375807019452 0.0808182811696644 0.078305419030476 0.0859551624513749 0.102173547873938 0.110373195222806 0.104887840961578 0.0866347967028697 0.0867226062085843 0.0924746447632758 0.102720828285174 0.0939522203893507 0.0823126744818889 0.0908582102478878 0.103376049094467 0.125023312465107 0.125906136298334 0.150623273183436 0.191912721361839 0.215469907246292 0.242460137370404 0.259146729169782 0.295578350419274 0.340937276992303 0.368209350081679 0.412603374240808 0.444591051989316 0.469596696872009 0.504061558216529 0.512590131445666 0.548443583014504 0.576591936267457 0.591751254212864 0.626972747123689 0.639506963823446 0.656717283460676 0.688818179484493 0.696304765770809 0.71240298688279 0.709769172003512 0.721147956359954 0.745943106193003 0.753551607682416 0.759405322707131 0.746609164293068 0.744720144805627 0.747880627791194 0.743553936113478 0.753191203632529 0.742692968929895 0.739550025554357 0.745808180721714 0.756791140248956 0.78072993441199 0.777068660791329 0.77946691140055 0.794180172202632 0.788505569984897 0.797583801688818 0.788845097523004 0.765716539002962 0.780326877748835 0.787010352811996 0.782372636180097 0.772159511336314 0.76937149195757 0.784535213524044 0.803755105405199 0.807614025468836 0.815335360826952 0.840251699198034 0.825837451541163 0.839866606140767 0.840570283808854 0.841490490525838 0.826607681535894 0.834919088156774 0.839590645684137 0.845193649196681 0.838505465644388 0.838724869188796 0.830752940992022 0.835392752447739 0.846319227488327 0.847792701014626 0.843677222317168 0.847466829779541 0.839343623742635 0.851013716812076 0.847911205203081 0.852387228173766 0.853570871216893 0.853684852360045 0.858322529478576 0.869297639794151 0.856461342536158 0.860841923489613 0.859709375168834 0.857456907165503 0.867530553822303 0.86537028202523 0.862916463489913 0.857125971550338 0.853082602155999 0.866759791074017 0.86159669799193 0.853139106847695 0.849981011684351 0.86196452415911 0.866492858067484 0.875032247895607 0.86503865854313 0.866779559762891 0.874848309230483 0.868473622169371 0.878199043104818 0.875077256431763 0.874087719727205 0.877288287099065 0.881171702581639 0.881265325494657 0.872392281190817 0.871997350084288 0.876162010046687 0.875467430630091 0.872284106065841 0.875323770651987 0.871498515306725 0.863490065750589 0.874388336720634 0.869917062484547 0.87414408225948 0.870200795132606 0.877007224933831 0.894633164347163 0.890074603288301 0.902993753486364 0.882293027051494 0.885350254344229 0.886560303478277 0.873960664285809 0.892271778062711 0.883659371029366 0.8880424842433 0.894853414946003 0.907039877886186 0.925333769530367 0.914534948193661 0.92167225011433 0.920023052490409 0.939932459508004 0.946407083021354 0.961108176644956 0.962631278702656 0.981698150548611 1];
rampDown=[0.948509532033036 0.988956079022078 0.992170230742876 0.988032223071807 0.97827288141617 0.985642254279523 1 0.99346191543377 0.967780850539307 0.96685802325231 0.960918106015062 0.965523781570725 0.964689897172557 0.961078185963817 0.961498242020093 0.964338236362123 0.948652129946432 0.953542443105194 0.965799240553875 0.961440462714719 0.96626525234358 0.966308877546851 0.947922359766621 0.944499227548375 0.94872281426011 0.961333257447823 0.945571691599879 0.957788316516863 0.964875836619889 0.915034358839554 0.94614354444659 0.954378537181412 0.93600322174368 0.947603867599536 0.945208161333888 0.931278256187503 0.923671326606471 0.915172790602474 0.913461521451027 0.908362518022036 0.915248379065542 0.873941950568248 0.846524651679065 0.815430332058461 0.792713296853549 0.764358737277334 0.739169026197824 0.684596869296654 0.639762711887849 0.613819690231547 0.547673719403686 0.510558942659553 0.486720903229241 0.4435477438309 0.41251072426503 0.361242789965736 0.341571297179933 0.319637507821993 0.297300930794428 0.282044249615756 0.274382033576828 0.26147893556643 0.250130472698177 0.235694964312543 0.207667696321296 0.191284622461747 0.171587546710908 0.160827699636723 0.162572883710862 0.16285774670382 0.188849197879136 0.213417132471127 0.20664702747559 0.183058997558463 0.185408194481197 0.190244891030617 0.178731029243586 0.16318466044596 0.151667877518727 0.139237410664801 0.12165994290534 0.120269118094493 0.119360507081589 0.1295362114234 0.1464188486931 0.155481954871526 0.168284350871153 0.150766037428772 0.147409545431902 0.156638837914038 0.15318418176615 0.139179092986554 0.121055511143841 0.106881789673515 0.0968973328498281 0.0903560118229136 0.0928235227446428 0.0976859105624813 0.105293917530223 0.122056634614371 0.13059606974835 0.150631607074657 0.137735933698567 0.123295357292588 0.118162644454614 0.10995417539872 0.0970046824315286 0.0846075451976215 0.0720898494658475 0.0704919262743557 0.0740294918209827 0.0822573185769725 0.0963633975114636 0.115614853858772 0.123209072197313 0.109871548604242 0.119279453768645 0.130942570306203 0.119418771072782 0.0997042004188851 0.0742829490374855 0.0647486970627031 0.0510925597365522 0.045863827088367 0.042952319872255 0.0495705449988617 0.0473105371585824 0.0518182085289661 0.060121547256883 0.0798596250993361 0.0897390072467838 0.0931905679266713 0.0877708926016275 0.0688881914183498 0.060659164341666 0.0505756925318924 0.0387385825912376 0.0338675266588315 0.0328687566924247 0.0323435691361318 0.0321076765635267 0.0352002250421378 0.0388429219811845 0.0433850862414412 0.0528673229263673 0.0539881134486781 0.0525543518586535 0.0432781621077257 0.033063972311327 0.030006748521193 0.0244963769951209 0.0199776712750445 0.0188545514081877 0.0171349220692097 0.0179341349472131 0.0167079117861384 0.0193932127965817 0.0237732948240338 0.0249987802886173 0.0215901186679346 0.0186780148036038 0.0165734233126092 0.0124896676765704 0.00823970837742318 0.0072081572691611 0.00840040332992698 0.00711550090457356 0.00338116175738272 0.00780656279877658 0.00600790209412821 0.00699427988053275 0.00724156801125731 0.00712159650666313 0.00865578514609755 0.00768633402657573 0.00603045211777831 0.00393923084037582 0.00113251973134941 1.98352603274928e-05 0];

inputs{1}= .2*stepInput
inputs{2}= .4*stepInput
inputs{3}= stepInput;
inputs{4}= rampUp;
inputs{5}= rampDown;
inputs{6}= complexRamp;
end

if nargin<6
    
    for j=1:numel(inputs)
        outputs{j}= repmat(mean(inputs{j}), 1, numel(inputs{j}));
    end
    
else
    

    figure;
 
    for j= 1:numel(inputs)
        
        subplot(2, numel(inputs), j);
       
        %func=rampSimHandle2(model, mediaInput,initialConditions, modelFeatures.options, outvar, gfpOutput);
        func=rampSimHandle2(model, inputs{j},initialConditions, modelFeatures.options, outvar, outputs{j});
        [x{j},y{j}]=func(params);
        plot(normalizeTS(inputs{j}), 'Color', [0 0 1]);
       subplot(2, numel(inputs), j+numel(inputs));
       numel([y{j}])
        plot(1:numel([y{j}]), normalizeTS([y{j}]') );
        
        
    end
        

    
    
    
    
    

lsqdiff=(vertcat(y{:}) -horzcat(outputs{:})').^2;
end








