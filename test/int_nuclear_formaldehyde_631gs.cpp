#include "../include/integral.h"

int test_int_nuclear_formaldehyde_631gs(int, char**) {
    // initialize the system
    System system("../example/molecule/formaldehyde.xyz", "6-31G*", 0, 1);

    // calculate the integral
    libint2::initialize();
    Matrix V = Integral::Nuclear(system);
    libint2::finalize();

    // create the expectation integral
    Matrix Vexp(system.shells.nbf(), system.shells.nbf()); Vexp << -64.61071832929831, -5.65517589772133, -0.09064364796294, -0.00229319450161, -0.01043675108548, -5.84190580994243, -0.01249053263203, -0.00031599810243, -0.00143816564023, -0.70510823216685, -0.00008664705890, -0.00038496460359, -0.70172415125713, -0.00000124769943, -0.70181104088592, -0.00003744068900, -0.68958432337121, 1.46254666979001, 0.03700097041781, 0.16839835221428, -1.72040572126900, 3.18689335293939, 0.08062521977389, 0.36694048838183, -1.98593566915684, -0.04986700545937, -0.22695378337404, -0.01609284039251, -0.00574142710416, -0.04096428666236, -0.00105212255578, -0.37995830806237, -0.00105212838453, -0.37995886133809, -5.65517589772133, -13.66067904943506, -0.61841001758019, -0.01564516084639, -0.07120401227395, -9.41937125330113, -0.40892452794348, -0.01034538548578, -0.04708375712839, -6.17846862491889, -0.00609388491990, -0.02707342870774, -5.94047093351971, -0.00008670893330, -5.94658709083037, -0.14473487656499, -2.11016479641508, 3.43368615389889, 0.08686883012348, 0.39535634010704, -3.49551162789840, 5.58615475526203, 0.14132413636978, 0.64319264811493, -4.49495846451955, -0.09972657356649, -0.45364838941999, -0.55649474977273, -0.01127235091083, -0.60728293935365, -0.02881031686374, -0.92574571671630, -0.02881041527748, -0.92574687848890, -0.09064364796294, -0.61841001758019, -14.94116053575950, -0.00682840255451, -0.03033759078996, -0.59995120655179, -5.52280141367957, -0.00691216680504, -0.03069748402368, -0.91207673151436, -0.01093045246125, -0.04917211285882, -0.24131896336447, -0.00028492807238, -0.24501280076209, -0.36431436345620, -2.44481346660678, 3.41421877843824, 0.11067563214871, 0.50401854435692, -1.83499963185664, 1.25319744217390, 0.07179097139476, 0.32711434207792, -4.66756888426559, -0.12733889498696, -0.57940952800175, -0.75694628971546, -0.01754423521872, -0.83650626046276, -0.04904270842088, -0.63044690979507, -0.05396793095627, -0.69181919713562, -0.00229319450161, -0.01564516084639, -0.00682840255451, -14.67447240634586, -0.00009800258554, -0.01517817122562, -0.00691216680504, -5.25289230804073, -0.00008779962897, -0.01093045246125, -0.24131896336447, -0.00028492807238, -0.01830917353206, -0.02777205522440, -0.00613872368958, -0.00921678966785, -0.06185134103581, 0.11067563573270, -0.95897406852773, 0.01303443302286, -0.04642366079055, 0.07179097571733, -1.58425208855617, 0.00861973803784, -0.17460892742544, 1.11431032736300, -0.01815996507476, 0.03743107244700, 0.12830539417381, -0.02121973972188, 0.00289804501799, 0.03562404414908, -0.00550412958242, -0.06907616682613, -0.01043675108548, -0.07120401227395, -0.03033759078996, -0.00009800258554, -14.68132155119667, -0.06907865626254, -0.03069748402368, -0.00008779962897, -5.25988253915101, -0.04917211285882, -0.00028492807238, -0.24501280076209, -0.02777205522440, -0.00613872368958, -0.08406922566106, -0.04194733713943, -0.28149703314477, 0.50401852582282, 0.01303443214130, -0.90523413580860, -0.21128275980468, 0.32711431968319, 0.00861973697358, -1.55021773720920, -0.79522458656496, -0.01817717548744, 1.03799571585125, -0.08716803766653, 0.02627418407939, 0.16149551733499, -0.02824060792979, -0.35413371742076, 0.01638000914059, 0.20188791635146, -5.84190580994243, -9.41937125330113, -0.59995120655179, -0.01517817122562, -0.06907865626254, -9.71900808660974, -1.06214468225704, -0.02687120934018, -0.12229581913063, -6.84890938840999, -0.01234309162239, -0.05477970966351, -6.36708349642952, -0.00012241297729, -6.37973797012060, -1.44419756019607, -3.91628433011064, 3.80466881058730, 0.09625431455930, 0.43807150048216, -5.24493255723057, 6.04066857795453, 0.15282287481853, 0.69552563445037, -4.71681252432856, -0.06476421745622, -0.29023902861384, -2.17709750845907, -0.00325631692020, -2.23041591911806, -0.38166026540360, -1.90851689363436, -0.38166090512503, -1.90851858807824, -0.01249053263203, -0.40892452794348, -5.52280141367957, -0.00691216680504, -0.03069748402368, -1.06214468225704, -8.05081061241346, -0.02505056337506, -0.11007992791002, -1.40859665200954, -0.01966985766667, -0.08796265373320, -0.31905620877490, -0.00060091142861, -0.32833286150767, -3.29373668662097, -5.83845897771221, 3.10904686498738, 0.14015475665636, 0.64289486324343, -4.80471399231863, 0.96711691758666, 0.11727803204725, 0.54057378270999, -4.68162651227037, -0.05687923674212, -0.24767594126506, -3.82326343490646, -0.00026310239685, -3.91950467785380, -0.89848141407129, -2.36097961782779, -0.99138671747658, -2.57120428273457, -0.00031599810243, -0.01034538548578, -0.00691216680504, -5.25289230804073, -0.00008779962897, -0.02687120934018, -0.02505056337506, -7.07745204139658, 0.00077181827575, -0.01966985766667, -0.31905620877490, -0.00060091142861, -0.02420045933336, -0.03669959256507, -0.00814408558058, -0.08332824968323, -0.14770718935056, 0.14015477441861, -2.44803529760481, 0.02081177249309, -0.12155447769468, 0.11727805281644, -3.69368449280600, 0.01984790620791, -0.18636833219755, 1.33921877531840, -0.00996044768504, -0.02884667304859, 0.15446176284359, -0.09910979150674, 0.05534124400136, 0.11693092078585, -0.10315322470075, -0.24171056490055, -0.00143816564023, -0.04708375712839, -0.03069748402368, -0.00008779962897, -5.25988253915101, -0.12229581913063, -0.11007992791002, 0.00077181827575, -7.10823958830387, -0.08796265373320, -0.00060091142861, -0.32833286150767, -0.03669959256507, -0.00814408558058, -0.11206488769869, -0.37924248181117, -0.67224304049401, 0.64289477239996, 0.02081176814968, -2.40152486526459, -0.55321710550786, 0.54057367559103, 0.01984790110653, -3.66694028198006, -0.84771958070572, -0.01023151367723, 1.29408983164731, -0.44020120145785, 0.03390550111147, -0.14263011091547, -0.52965048466261, -1.23624826042611, 0.31205120144851, 0.66835578576803, -0.70510823216685, -6.17846862491889, -0.91207673151436, -0.01093045246125, -0.04917211285882, -6.84890938840999, -1.40859665200954, -0.01966985766667, -0.08796265373320, -9.90672673556426, -0.01651196918033, -0.07318164000658, -3.08765594033528, -0.00060447169786, -3.09591231404099, -2.04704181163176, -4.63843503101378, 5.22863313173207, 0.20512698989977, 0.93582735526197, -4.21927420938442, 4.14023531177360, 0.15107823117773, 0.68940673701340, -5.63662683063952, -0.15455067208536, -0.70249448209090, -2.09698556223253, -0.02543367759738, -2.22534142832738, -0.20314558428077, -1.43900703201508, -0.23873068553741, -1.54566733778784, -0.00008664705890, -0.00609388491990, -0.01093045246125, -0.24131896336447, -0.00028492807238, -0.01234309162239, -0.01966985766667, -0.31905620877490, -0.00060091142861, -0.01651196918033, -3.08765594033528, -0.00060447169786, -0.00885421946085, -0.01303538688779, -0.00303260418244, -0.05094656147316, -0.08375480654180, 0.12062168943646, -1.44040971043774, 0.01916025359670, -0.03812671359542, 0.02236382215564, -0.91833533054671, 0.00611276289698, -0.15226186076358, 1.27602145916043, -0.01779123339180, 0.03896697085444, 0.18833169098623, -0.04676312438145, 0.00966940303548, 0.02545857607897, -0.02113441592232, -0.06687084852238, -0.00038496460359, -0.02707342870774, -0.04917211285882, -0.00028492807238, -0.24501280076209, -0.05477970966351, -0.08796265373320, -0.00060091142861, -0.32833286150767, -0.07318164000658, -0.00060447169786, -3.09591231404099, -0.01303538688779, -0.00303260418244, -0.03991482992108, -0.23185333358445, -0.38049815293790, 0.54948006592779, 0.01911603532721, -1.36727922108721, -0.17244737228475, 0.10158783499787, 0.00606177277371, -0.89988952253988, -0.69352912018455, -0.01782795356304, 1.20403089119642, -0.19984560054443, 0.03911661632975, 0.16678898142107, -0.09815724303603, -0.31134219179119, 0.06099059313645, 0.16567773126835, -0.70172415125713, -5.94047093351971, -0.24131896336447, -0.01830917353206, -0.02777205522440, -6.36708349642952, -0.31905620877490, -0.02420045933336, -0.03669959256507, -3.08765594033528, -0.00885421946085, -0.01303538688779, -8.91774872783473, 0.00009445478075, -2.97591921972013, -0.03461089171656, -1.33278097956570, 1.90705204792055, -0.02483637559366, 0.21963202536751, -2.71762033818449, 4.17769356056259, 0.05916725441143, 0.48106459458277, -2.18141940894593, 0.03682941915047, -0.20957310303046, -1.07686587624370, 0.00476901486109, -0.38576009868581, -0.02464455423156, -0.76027280062317, -0.02618047201786, -0.76487714721194, -0.00000124769943, -0.00008670893330, -0.00028492807238, -0.02777205522440, -0.00613872368958, -0.00012241297729, -0.00060091142861, -0.03669959256507, -0.00814408558058, -0.00060447169786, -0.01303538688779, -0.00303260418244, 0.00009445478075, -2.97591921972013, 0.00013947213241, -0.00585297245696, -0.00900567513936, 0.01753421892877, -0.16583752612911, -0.03453784398347, -0.00339033137793, 0.00423767169799, -0.10574214853395, -0.02278653450975, -0.02652056948818, 0.18829751740944, 0.03896535828469, 0.00501426907177, -0.33777629746702, 0.00481525206882, 0.00603176949120, 0.01665619907729, 0.00661497401201, 0.01840421870192, -0.70181104088592, -5.94658709083037, -0.24501280076209, -0.00613872368958, -0.08406922566106, -6.37973797012060, -0.32833286150767, -0.00814408558058, -0.11206488769869, -3.09591231404099, -0.00303260418244, -0.03991482992108, -2.97591921972013, 0.00013947213241, -8.93816274977703, -0.06008456441472, -1.37774356896683, 1.98851989862848, 0.05054265419467, -0.10489381257754, -2.74163680448330, 4.20575375231485, 0.10659117039763, 0.27151177564075, -2.30187667755509, -0.04872973051856, 0.15785778497076, -0.38468971429020, 0.00460671764286, -1.04600115409973, -0.08177420586064, -0.92817695803401, -0.04465397454275, -0.81691552660220, -0.00003744068900, -0.14473487656499, -0.36431436345620, -0.00921678966785, -0.04194733713943, -1.44419756019607, -3.29373668662097, -0.08332824968323, -0.37924248181117, -2.04704181163176, -0.05094656147316, -0.23185333358445, -0.03461089171656, -0.00585297245696, -0.06008456441472, -38.59305990970996, -3.38555861569718, 0.11900328172155, 0.00301066150394, 0.01370211025570, -4.07106242610478, 0.01999183171669, 0.00050577292720, 0.00230187166633, -1.18620005381586, -0.00054107561607, -0.00167712944419, -1.16806168675990, 0.00066841456213, -1.17198778722858, -0.63612646359122, -1.98521538170805, -0.63612644738161, -1.98521536916219, -0.68958432337121, -2.11016479641508, -2.44481346660678, -0.06185134103581, -0.28149703314477, -3.91628433011064, -5.83845897771221, -0.14770718935056, -0.67224304049401, -4.63843503101378, -0.08375480654180, -0.38049815293790, -1.33278097956570, -0.00900567513936, -1.37774356896683, -3.38555861569718, -10.15593952482487, 0.85552014342860, 0.02164378624859, 0.09850510855273, -7.71265971305316, 0.59954177285671, 0.01516779535720, 0.06903160091189, -7.05828983207658, -0.01708684960682, -0.05348426240199, -6.48334275236839, 0.02062295249507, -6.60489898507538, -2.22192643027363, -4.20376785257915, -2.22192656789909, -4.20376805101338, 1.46254666979001, 3.43368615389889, 3.41421877843824, 0.11067563573270, 0.50401852582282, 3.80466881058730, 3.10904686498738, 0.14015477441861, 0.64289477239996, 5.22863313173207, 0.12062168943646, 0.54948006592779, 1.90705204792055, 0.01753421892877, 1.98851989862848, 0.11900328172155, 0.85552014342860, -11.07988477924484, -0.02020119541587, -0.06354555127049, 0.83447283589656, -5.23850610521795, -0.01860285709485, -0.06029614785319, 1.40826200776045, 0.01658039551991, 0.09905535305982, 0.32802590566612, 0.01208687046013, 0.27934234862894, -1.23009994715514, -0.59118032676472, -1.75835592922333, -0.97443762343393, 0.03700097041781, 0.08686883012348, 0.11067563214871, -0.95897406852773, 0.01303443214130, 0.09625431455930, 0.14015475665636, -2.44803529760481, 0.02081176814968, 0.20512698989977, -1.44040971043774, 0.01911603532721, -0.02483637559366, -0.16583752612911, 0.05054265419467, 0.00301066150394, 0.02164378624859, -0.02020119541587, -10.39885554092802, 0.02409065580256, 0.02111131147240, -0.01860285709485, -4.60403521611573, 0.02053003635102, 0.01658039551991, 0.32802590566612, 0.01208687046013, 0.02488712851637, 0.03832446203515, 0.00952586400963, 0.41280158243891, 0.30711554026559, -0.48840692659506, -0.34672440537368, 0.16839835221428, 0.39535634010704, 0.50401854435692, 0.01303443302286, -0.90523413580860, 0.43807150048216, 0.64289486324343, 0.02081177249309, -2.40152486526459, 0.93582735526197, 0.01916025359670, -1.36727922108721, 0.21963202536751, -0.03453784398347, -0.10489381257754, 0.01370211025570, 0.09850510855273, -0.06354555127049, 0.02409065580256, -10.54111046890088, 0.09608170582815, -0.06029614785319, 0.02053003635102, -4.72675601050454, 0.09905535305982, 0.01208687046013, 0.27934234862894, 0.03832446203515, 0.00952586400963, 0.09470107602406, -2.56502798501253, -1.82627631355546, 2.22093762018872, 1.64601172145414, -1.72040572126900, -3.49551162789840, -1.83499963185664, -0.04642366079055, -0.21128275980468, -5.24493255723057, -4.80471399231863, -0.12155447769468, -0.55321710550786, -4.21927420938442, -0.03812671359542, -0.17244737228475, -2.71762033818449, -0.00339033137793, -2.74163680448330, -4.07106242610478, -7.71265971305316, 0.83447283589656, 0.02111131147240, 0.09608170582815, -8.08419705680052, 1.07845012326693, 0.02728369031712, 0.12417338654602, -5.95562192490182, -0.01890801141763, -0.06017093090357, -5.31533330264414, 0.02190350094884, -5.44525445326547, -2.65841424093863, -4.94772563851828, -2.65841446301488, -4.94772602615583, 3.18689335293939, 5.58615475526203, 1.25319744217390, 0.07179097571733, 0.32711431968319, 6.04066857795453, 0.96711691758666, 0.11727805281644, 0.54057367559103, 4.14023531177360, 0.02236382215564, 0.10158783499787, 4.17769356056259, 0.00423767169799, 4.20575375231485, 0.01999183171669, 0.59954177285671, -5.23850610521795, -0.01860285709485, -0.06029614785319, 1.07845012326693, -7.10691631300662, -0.03430330625044, -0.11677879041180, 1.17829728463029, 0.01526962414614, 0.09246086738291, 0.24026828468116, 0.01189491380525, 0.19344634513709, -1.42641339192879, -1.01017100772038, -2.03175157527765, -1.60334672659583, 0.08062521977389, 0.14132413636978, 0.07179097139476, -1.58425208855617, 0.00861973697358, 0.15282287481853, 0.11727803204725, -3.69368449280600, 0.01984790110653, 0.15107823117773, -0.91833533054671, 0.00606177277371, 0.05916725441143, -0.10574214853395, 0.10659117039763, 0.00050577292720, 0.01516779535720, -0.01860285709485, -4.60403521611573, 0.02053003635102, 0.02728369031712, -0.03430330625044, -5.91392260456540, 0.03265265439844, 0.01526962414614, 0.24026828468116, 0.01189491380525, 0.01822538731079, 0.02820519315322, 0.00728724583999, 0.47261130636860, 0.47292128321794, -0.56009989427178, -0.53904113627911, 0.36694048838183, 0.64319264811493, 0.32711434207792, 0.00861973803784, -1.55021773720920, 0.69552563445037, 0.54057378270999, 0.01984790620791, -3.66694028198006, 0.68940673701340, 0.00611276289698, -0.89988952253988, 0.48106459458277, -0.02278653450975, 0.27151177564075, 0.00230187166633, 0.06903160091189, -0.06029614785319, 0.02053003635102, -4.72675601050454, 0.12417338654602, -0.11677879041180, 0.03265265439844, -6.11417660701308, 0.09246086738291, 0.01189491380525, 0.19344634513709, 0.02820519315322, 0.00728724583999, 0.06494196235822, -2.94124872282969, -2.83752789012994, 2.54307618490103, 2.53660843313942, -1.98593566915684, -4.49495846451955, -4.66756888426559, -0.17460892742544, -0.79522458656496, -4.71681252432856, -4.68162651227037, -0.18636833219755, -0.84771958070572, -5.63662683063952, -0.15226186076358, -0.69352912018455, -2.18141940894593, -0.02652056948818, -2.30187667755509, -1.18620005381586, -7.05828983207658, 1.40826200776045, 0.01658039551991, 0.09905535305982, -5.95562192490182, 1.17829728463029, 0.01526962414614, 0.09246086738291, -9.85467677039385, -0.02587476044900, -0.08080289736703, -2.99519156814214, 0.00981040187882, -3.05922748526586, -1.40928016378648, -3.11661651711112, -1.89939266860629, -3.21162679662954, -0.04986700545937, -0.09972657356649, -0.12733889498696, 1.11431032736300, -0.01817717548744, -0.06476421745622, -0.05687923674212, 1.33921877531840, -0.01023151367723, -0.15455067208536, 1.27602145916043, -0.01782795356304, 0.03682941915047, 0.18829751740944, -0.04872973051856, -0.00054107561607, -0.01708684960682, 0.01658039551991, 0.32802590566612, 0.01208687046013, -0.01890801141763, 0.01526962414614, 0.24026828468116, 0.01189491380525, -0.02587476044900, -2.99519156814214, 0.00981040187882, -0.01384373956603, -0.01422513209596, -0.00560395458554, 0.15698215114194, 0.02240942280664, -0.26728525630397, -0.05983638409061, -0.22695378337404, -0.45364838941999, -0.57940952800175, -0.01815996507476, 1.03799571585125, -0.29023902861384, -0.24767594126506, -0.00996044768504, 1.29408983164731, -0.70249448209090, -0.01779123339180, 1.20403089119642, -0.20957310303046, 0.03896535828469, 0.15785778497076, -0.00167712944419, -0.05348426240199, 0.09905535305982, 0.01208687046013, 0.27934234862894, -0.06017093090357, 0.09246086738291, 0.01189491380525, 0.19344634513709, -0.08080289736703, 0.00981040187882, -3.05922748526586, -0.01422513209596, -0.00560395458554, -0.04696263481275, -0.98279631103660, -0.19674548022015, 1.20918183449747, 0.22817757079262, -0.01609284039251, -0.55649474977273, -0.75694628971546, 0.03743107244700, -0.08716803766653, -2.17709750845907, -3.82326343490646, -0.02884667304859, -0.44020120145785, -2.09698556223253, 0.03896697085444, -0.19984560054443, -1.07686587624370, 0.00501426907177, -0.38468971429020, -1.16806168675990, -6.48334275236839, 0.32802590566612, 0.02488712851637, 0.03832446203515, -5.31533330264414, 0.24026828468116, 0.01822538731079, 0.02820519315322, -2.99519156814214, -0.01384373956603, -0.01422513209596, -8.51623530180116, 0.01608742194507, -2.87212105551323, -0.96533557331072, -2.83840677291140, -0.98648911272752, -2.84250753899627, -0.00574142710416, -0.01127235091083, -0.01754423521872, 0.12830539417381, 0.02627418407939, -0.00325631692020, -0.00026310239685, 0.15446176284359, 0.03390550111147, -0.02543367759738, 0.18833169098623, 0.03911661632975, 0.00476901486109, -0.33777629746702, 0.00460671764286, 0.00066841456213, 0.02062295249507, 0.01208687046013, 0.03832446203515, 0.00952586400963, 0.02190350094884, 0.01189491380525, 0.02820519315322, 0.00728724583999, 0.00981040187882, -0.01422513209596, -0.00560395458554, 0.01608742194507, -2.87212105551323, 0.02868780521215, 0.32846988333996, 0.09092614074035, 0.33650259071400, 0.09248334373609, -0.04096428666236, -0.60728293935365, -0.83650626046276, -0.02121973972188, 0.16149551733499, -2.23041591911806, -3.91950467785380, -0.09910979150674, -0.14263011091547, -2.22534142832738, -0.04676312438145, 0.16678898142107, -0.38576009868581, 0.00481525206882, -1.04600115409973, -1.17198778722858, -6.60489898507538, 0.27934234862894, 0.00952586400963, 0.09470107602406, -5.44525445326547, 0.19344634513709, 0.00728724583999, 0.06494196235822, -3.05922748526586, -0.00560395458554, -0.04696263481275, -2.87212105551323, 0.02868780521215, -8.77404265644387, -2.95447556903124, -3.36898669488541, -2.44320990750352, -3.26987619758297, -0.00105212255578, -0.02881031686374, -0.04904270842088, 0.00289804501799, -0.02824060792979, -0.38166026540360, -0.89848141407129, 0.05534124400136, -0.52965048466261, -0.20314558428077, 0.00966940303548, -0.09815724303603, -0.02464455423156, 0.00603176949120, -0.08177420586064, -0.63612646359122, -2.22192643027363, -1.23009994715514, 0.41280158243891, -2.56502798501253, -2.65841424093863, -1.42641339192879, 0.47261130636860, -2.94124872282969, -1.40928016378648, 0.15698215114194, -0.98279631103660, -0.96533557331072, 0.32846988333996, -2.95447556903124, -6.89987672744151, -4.18304908550170, -0.11365574859280, -0.95070802886904, -0.37995830806237, -0.92574571671630, -0.63044690979507, 0.03562404414908, -0.35413371742076, -1.90851689363436, -2.36097961782779, 0.11693092078585, -1.23624826042611, -1.43900703201508, 0.02545857607897, -0.31134219179119, -0.76027280062317, 0.01665619907729, -0.92817695803401, -1.98521538170805, -4.20376785257915, -0.59118032676472, 0.30711554026559, -1.82627631355546, -4.94772563851828, -1.01017100772038, 0.47292128321794, -2.83752789012994, -3.11661651711112, 0.02240942280664, -0.19674548022015, -2.83840677291140, 0.09092614074035, -3.36898669488541, -4.18304908550170, -5.61089940365182, -0.95070808991821, -2.47876861470137, -0.00105212838453, -0.02881041527748, -0.05396793095627, -0.00550412958242, 0.01638000914059, -0.38166090512503, -0.99138671747658, -0.10315322470075, 0.31205120144851, -0.23873068553741, -0.02113441592232, 0.06099059313645, -0.02618047201786, 0.00661497401201, -0.04465397454275, -0.63612644738161, -2.22192656789909, -1.75835592922333, -0.48840692659506, 2.22093762018872, -2.65841446301488, -2.03175157527765, -0.56009989427178, 2.54307618490103, -1.89939266860629, -0.26728525630397, 1.20918183449747, -0.98648911272752, 0.33650259071400, -2.44320990750352, -0.11365574859280, -0.95070808991821, -6.89987734736981, -4.18304949445143, -0.37995886133809, -0.92574687848890, -0.69181919713562, -0.06907616682613, 0.20188791635146, -1.90851858807824, -2.57120428273457, -0.24171056490055, 0.66835578576803, -1.54566733778784, -0.06687084852238, 0.16567773126835, -0.76487714721194, 0.01840421870192, -0.81691552660220, -1.98521536916219, -4.20376805101338, -0.97443762343393, -0.34672440537368, 1.64601172145414, -4.94772602615583, -1.60334672659583, -0.53904113627911, 2.53660843313942, -3.21162679662954, -0.05983638409061, 0.22817757079262, -2.84250753899627, 0.09248334373609, -3.26987619758297, -0.95070802886904, -2.47876861470137, -4.18304949445143, -5.61090001835214;

    // print the results
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR: " << V << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "COMPUTED NUCLEAR NORM: " << V.norm() << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR: " << Vexp << std::endl;
    std::cout << std::fixed << std::setprecision(14) << "EXPECTED NUCLEAR NORM: " << Vexp.norm() << std::endl;

    // return success or failure based on the error
    return (V - Vexp).norm() > 1e-8;
}
