#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Archimedes algorithm in the formulation of Pfaff.

Description:
The programmed iterative Archimedes algorithm uses an approach developed
by the author with the intention to improve the calculation result of Pi.

The algorithm by Pfaff is used and approaches by Snellius and Dörrie are
considered in the development. In general it is a mixture of harmonic,
geometric and arithmetic mean in conjunction with a Snellius acceleration
in squared form. The acceleration of the convergence of the calculation
is in principle a suitable selected weighted arithmetic mean.

In the presented form, the precision and the number of iterations can be
predicted via the number of desired decimal places of Pi.

The algorithm is tested up to 70.000 correct places of Pi in a manageable
period of time so far. Increasing the precision slows down the calculation
process.

One advantage of the Pfaff algorithm is that we do not have to worry about
whether the calculation exceeds the allowed range of values.

Calculation results:
As it seems the radius r has no influence on precision and iterations.

Ludolph van Ceulen calculated 35 places of Pi applying the approach of
Archimedes using 6*2^62 edges. We can do so with only 6*2^18 edges.
This is significant faster.

Preliminary conclusion:
Using the weighted arithmetic mean as done in NETZ gives only good
results in combination with DÖRRIE.

Spin-off:
Changing the formula NETZ for calculating the Archimedes constants
resulted in:

Weight | Name      | Precision | Iterations | Correct Places
------------------------------------------------------------
 2/3   | Snellius  | 1002      | 828        | 1000
 3/4   |           | 1002      | 828        | 1000
 4/5   | Netz      | 1002      | 552        | 1000            <-- beacon
 5/6   |           | 1002      | 828        | 1000
 6/7   |           | 1002      | 829        | 1000
 etc.  |           | 1002      | 829        | 1000

2/3 means a weight of 2 for one element and therefor the sum is of n
elements to 3.

One sees that the required precision stays the same for a requested number
of places. The number of required iterations changes only for a weight of
4/5 in the arithmetic mean. This is an astonishing result.

Bugs:
No bugs known yet.

Limitations:
No limitations known yet.

To-Do:
1. Check if calculations can be done on high-end graphic cards. So far the
standard Python module decimal stands against that.
2. Check the influence of rounding using the decimal module.
3. Improve code and documentation.

Test system:
No high-end standard Linux PC: Python 3.8.10; Linux Mint 20.3 Una, Ubuntu
Focal, GNU/Linux, x86_64

Acknowledgement:
Eve Andersson
www.eveandersson.com/pi/digits/
'''
# pylint: disable=invalid-name
# pylint: disable=global-statement
# pylint: disable=protected-access

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.3"

# Import some standard Python modules.
import sys
import os

# Import names from the standard Python module decimal.
from decimal import Decimal as D
from decimal import getcontext, ROUND_HALF_DOWN, ROUND_UP

# Set some user defined constants.
RADIUS = 1        # radius of the circle
PLACES = 10000    # number of requested places
PROGRESS = True   # show or hide calculation progress

# Choose the calculation method:
# WEIGHTED-ARITHMETIC-MEAN -> 4
# ARITHMETIC-MEAN          -> 3
# SNELLIUS                 -> 2
# DOERRIE                  -> 1
# NETZ                     -> 0
# Warning: Set OVERRUN to True and use userdefinded precision and
#          iteration. Precalculated values are only valid for NETZ.
METHOD = 0

# Overrun the calculation of precision and iteration.
OVERRUN = False

# Check value of constant OVERRUN. Set the constants.
# Ludolph van Ceulen:
# PLACES, PRECISION, ITERATION = 35, 37, 18
# Dörrie:
# PLACES, PRECISION, ITERATION = 1000, 1000, 829
# Snellius:
# PLACES, PRECISION, ITERATION = 1000, 1002, 830
# Arithmetic mean:
# PLACES, PRECISION, ITERATION = 1000, 1003, 1660
# Weighted arithmetic mean:
# PLACES, PRECISION, ITERATION = 1000, 1002, 1659
if OVERRUN:
    # Set the user defined constants.
    PLACES, PRECISION, ITERATION = 100, 102, 54

# Calculate the initial values.
def init_values(places, offset=16):
    '''Predict precision and iteration by places.

    Base values developed from data observations.
    '''
    # Set the base precision and base iteration.
    baseprec, baseiter = 1.02, 0.56
    # Calculate precision and iteration.
    calcprec = baseprec*places
    calciter = baseiter*places
    # Round them up.
    precision = D(calcprec).quantize(D('1'), rounding=ROUND_UP)
    iteration = D(calciter).quantize(D('1'), rounding=ROUND_UP)
    # Consider the offset.
    precision += D(offset)
    iteration += D(offset)/D(2)
    # Return required precision and iteration.
    return int(precision), int(iteration)

# Check value of constant OVERRUN.
if not OVERRUN:
    # Initialise the script defined constants.
    # Offset is a kind of minimum precision.
    PRECISION, ITERATION = init_values(PLACES, offset=4)

# Set the precision and the rounding method.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_DOWN

# Define a heredoc consisting of Pi with 10000 places.
PI_HEREDOC = '''
3.
1415926535897932384626433832795028841971693993751058209749445923078164062862089
9862803482534211706798214808651328230664709384460955058223172535940812848111745
0284102701938521105559644622948954930381964428810975665933446128475648233786783
1652712019091456485669234603486104543266482133936072602491412737245870066063155
8817488152092096282925409171536436789259036001133053054882046652138414695194151
1609433057270365759591953092186117381932611793105118548074462379962749567351885
7527248912279381830119491298336733624406566430860213949463952247371907021798609
4370277053921717629317675238467481846766940513200056812714526356082778577134275
7789609173637178721468440901224953430146549585371050792279689258923542019956112
1290219608640344181598136297747713099605187072113499999983729780499510597317328
1609631859502445945534690830264252230825334468503526193118817101000313783875288
6587533208381420617177669147303598253490428755468731159562863882353787593751957
7818577805321712268066130019278766111959092164201989380952572010654858632788659
3615338182796823030195203530185296899577362259941389124972177528347913151557485
7242454150695950829533116861727855889075098381754637464939319255060400927701671
1390098488240128583616035637076601047101819429555961989467678374494482553797747
2684710404753464620804668425906949129331367702898915210475216205696602405803815
0193511253382430035587640247496473263914199272604269922796782354781636009341721
6412199245863150302861829745557067498385054945885869269956909272107975093029553
2116534498720275596023648066549911988183479775356636980742654252786255181841757
4672890977772793800081647060016145249192173217214772350141441973568548161361157
3525521334757418494684385233239073941433345477624168625189835694855620992192221
8427255025425688767179049460165346680498862723279178608578438382796797668145410
0953883786360950680064225125205117392984896084128488626945604241965285022210661
1863067442786220391949450471237137869609563643719172874677646575739624138908658
3264599581339047802759009946576407895126946839835259570982582262052248940772671
9478268482601476990902640136394437455305068203496252451749399651431429809190659
2509372216964615157098583874105978859597729754989301617539284681382686838689427
7415599185592524595395943104997252468084598727364469584865383673622262609912460
8051243884390451244136549762780797715691435997700129616089441694868555848406353
4220722258284886481584560285060168427394522674676788952521385225499546667278239
8645659611635488623057745649803559363456817432411251507606947945109659609402522
8879710893145669136867228748940560101503308617928680920874760917824938589009714
9096759852613655497818931297848216829989487226588048575640142704775551323796414
5152374623436454285844479526586782105114135473573952311342716610213596953623144
2952484937187110145765403590279934403742007310578539062198387447808478489683321
4457138687519435064302184531910484810053706146806749192781911979399520614196634
2875444064374512371819217999839101591956181467514269123974894090718649423196156
7945208095146550225231603881930142093762137855956638937787083039069792077346722
1825625996615014215030680384477345492026054146659252014974428507325186660021324
3408819071048633173464965145390579626856100550810665879699816357473638405257145
9102897064140110971206280439039759515677157700420337869936007230558763176359421
8731251471205329281918261861258673215791984148488291644706095752706957220917567
1167229109816909152801735067127485832228718352093539657251210835791513698820914
4421006751033467110314126711136990865851639831501970165151168517143765761835155
6508849099898599823873455283316355076479185358932261854896321329330898570642046
7525907091548141654985946163718027098199430992448895757128289059232332609729971
2084433573265489382391193259746366730583604142813883032038249037589852437441702
9132765618093773444030707469211201913020330380197621101100449293215160842444859
6376698389522868478312355265821314495768572624334418930396864262434107732269780
2807318915441101044682325271620105265227211166039666557309254711055785376346682
0653109896526918620564769312570586356620185581007293606598764861179104533488503
4611365768675324944166803962657978771855608455296541266540853061434443185867697
5145661406800700237877659134401712749470420562230538994561314071127000407854733
2699390814546646458807972708266830634328587856983052358089330657574067954571637
7525420211495576158140025012622859413021647155097925923099079654737612551765675
1357517829666454779174501129961489030463994713296210734043751895735961458901938
9713111790429782856475032031986915140287080859904801094121472213179476477726224
1425485454033215718530614228813758504306332175182979866223717215916077166925474
8738986654949450114654062843366393790039769265672146385306736096571209180763832
7166416274888800786925602902284721040317211860820419000422966171196377921337575
1149595015660496318629472654736425230817703675159067350235072835405670403867435
1362222477158915049530984448933309634087807693259939780541934144737744184263129
8608099888687413260472156951623965864573021631598193195167353812974167729478672
4229246543668009806769282382806899640048243540370141631496589794092432378969070
6977942236250822168895738379862300159377647165122893578601588161755782973523344
6042815126272037343146531977774160319906655418763979293344195215413418994854447
3456738316249934191318148092777710386387734317720754565453220777092120190516609
6280490926360197598828161332316663652861932668633606273567630354477628035045077
7235547105859548702790814356240145171806246436267945612753181340783303362542327
8394497538243720583531147711992606381334677687969597030983391307710987040859133
7464144282277263465947047458784778720192771528073176790770715721344473060570073
3492436931138350493163128404251219256517980694113528013147013047816437885185290
9285452011658393419656213491434159562586586557055269049652098580338507224264829
3972858478316305777756068887644624824685792603953527734803048029005876075825104
7470916439613626760449256274204208320856611906254543372131535958450687724602901
6187667952406163425225771954291629919306455377991403734043287526288896399587947
5729174642635745525407909145135711136941091193932519107602082520261879853188770
5842972591677813149699009019211697173727847684726860849003377024242916513005005
1683233643503895170298939223345172201381280696501178440874519601212285993716231
3017114448464090389064495444006198690754851602632750529834918740786680881833851
0228334508504860825039302133219715518430635455007668282949304137765527939751754
6139539846833936383047461199665385815384205685338621867252334028308711232827892
1250771262946322956398989893582116745627010218356462201349671518819097303811980
0497340723961036854066431939509790190699639552453005450580685501956730229219139
3391856803449039820595510022635353619204199474553859381023439554495977837790237
4216172711172364343543947822181852862408514006660443325888569867054315470696574
7458550332323342107301545940516553790686627333799585115625784322988273723198987
5714159578111963583300594087306812160287649628674460477464915995054973742562690
1049037781986835938146574126804925648798556145372347867330390468838343634655379
4986419270563872931748723320837601123029911367938627089438799362016295154133714
2489283072201269014754668476535761647737946752004907571555278196536213239264061
6013635815590742202020318727760527721900556148425551879253034351398442532234157
6233610642506390497500865627109535919465897514131034822769306247435363256916078
1547818115284366795706110861533150445212747392454494542368288606134084148637767
0096120715124914043027253860764823634143346235189757664521641376796903149501910
8575984423919862916421939949072362346468441173940326591840443780513338945257423
9950829659122850855582157250310712570126683024029295252201187267675622041542051
6184163484756516999811614101002996078386909291603028840026910414079288621507842
4516709087000699282120660418371806535567252532567532861291042487761825829765157
9598470356222629348600341587229805349896502262917487882027342092222453398562647
6691490556284250391275771028402799806636582548892648802545661017296702664076559
0429099456815065265305371829412703369313785178609040708667114965583434347693385
7817113864558736781230145876871266034891390956200993936103102916161528813843790
9904231747336394804575931493140529763475748119356709110137751721008031559024853
0906692037671922033229094334676851422144773793937517034436619910403375111735471
9185504644902636551281622882446257591633303910722538374218214088350865739177150
9682887478265699599574490661758344137522397096834080053559849175417381883999446
9748676265516582765848358845314277568790029095170283529716344562129640435231176
0066510124120065975585127617858382920419748442360800719304576189323492292796501
9875187212726750798125547095890455635792122103334669749923563025494780249011419
5212382815309114079073860251522742995818072471625916685451333123948049470791191
5326734302824418604142636395480004480026704962482017928964766975831832713142517
0296923488962766844032326092752496035799646925650493681836090032380929345958897
0695365349406034021665443755890045632882250545255640564482465151875471196218443
9658253375438856909411303150952617937800297412076651479394259029896959469955657
6121865619673378623625612521632086286922210327488921865436480229678070576561514
4632046927906821207388377814233562823608963208068222468012248261177185896381409
1839036736722208883215137556003727983940041529700287830766709444745601345564172
5437090697939612257142989467154357846878861444581231459357198492252847160504922
1242470141214780573455105008019086996033027634787081081754501193071412233908663
9383395294257869050764310063835198343893415961318543475464955697810382930971646
5143840700707360411237359984345225161050702705623526601276484830840761183013052
7932054274628654036036745328651057065874882256981579367897669742205750596834408
6973502014102067235850200724522563265134105592401902742162484391403599895353945
9094407046912091409387001264560016237428802109276457931065792295524988727584610
1264836999892256959688159205600101655256375678
'''

# ----------------------------------------------------------------------
# Helper function remove_whitestrings()
# ----------------------------------------------------------------------
def remove_whitespaces(string):
    '''Remove all whitespaces defined by a list from a string.'''
    # Define the list with whitespaces to remove.
    mapping = [("\n", ""), ("\r", ""), ("\t", ""), (" ", "")]
    # Remove the whitespaces from the given string.
    for k, v in mapping:
        string = string.replace(k, v)
    # Return the trimmed string.
    return string

# *********************************
# Generator function chunk_string()
# *********************************
def chunk_string(string, length):
    '''Generate and return a chunk string.'''
    # Return the chunk string.
    return (string[0+i:length+i] for i in range(0, len(string), length))

# ----------------------------------------------------------------------
# Helper function print_pi()
# ----------------------------------------------------------------------
def print_pi(prtnum, prtlen, maxln=-1):
    '''Print the circle number Pi in chunks.'''
    # Initialise the local variables.
    outstr = str(prtnum)
    outlen = int(prtlen)
    count = 1
    # Print the leading 3 with decimal point to first line.
    print("{:<4d}".format(0), "3.")
    # Run over the places of Pi from the given string.
    for i in chunk_string(outstr[2:], outlen):
        # Print a line until there is data.
        if maxln != -1 and maxln < count:
            break
        # Print a line with the places.
        print("{:<4d}".format(count), i)
        # Increment the line counter.
        count += 1
    # Print an empty line.
    print("\r")
    # End of function. Return 1.
    return 1

# ----------------------------------------------------------------------
# Helper function correct_digits()
# ----------------------------------------------------------------------
def correct_digits(chkpi, refpi):
    '''Calculate the correct digits of a given pi number.'''
    # Initialise the local variables.
    print(len(chkpi), len(refpi))
    correct = ''
    idx = 0
    # Run over the digits of Pi.
    for char in str(refpi):
        # Exit condition.
        try:
            # Compare the calculated value with the reference value.
            if char == str(chkpi)[idx]:
                # Add the correct char to string.
                correct += char
                # Increment counter.
                idx += 1
            else:
                # Leave loop.
                break
        except:
            pass    # Return the correct digits and thenumber of correct digits.
    return (correct, idx-2)

# ----------------------------------------------------------------------
# Function hide_cursor()
# ----------------------------------------------------------------------
def hide_cursor():
    '''Hide the cursor.'''
    sys.stdout.write("\x1b[?25l")
    sys.stdout.flush()

# ----------------------------------------------------------------------
# Function show_cursor()
# ----------------------------------------------------------------------
def show_cursor():
    '''Show the cursor.'''
    sys.stdout.write("\x1b[?25h")
    sys.stdout.flush()

# ----------------------------------------------------------------------
# Function print_iteration()
# ----------------------------------------------------------------------
def print_iteration(citer):
    '''Print progress in form of the current iteration.'''
    if citer % 100 == 0 and citer >= 100:
        string = "Iteration: " + str(citer)
        sys.stdout.write(string)
        sys.stdout.write("\r")
        sys.stdout.flush()

# ----------------------------------------------------------------------
# Function archimedes_weighted_arithmetic_mean()
# ----------------------------------------------------------------------
def archimedes_weighted_arithmetic_mean(a1, b1, r):
    '''Function Archimedes-Arithmetic-Mean.

         a₁ + w⋅b₁
    ac = ─────────
         (w + 1)⋅r
    '''
    # Print method which is used.
    print("*** WEIGHTED-ARITHMETIC-MEAN ***")
    # Calculate the Archimedes constant.
    w = 4
    ac = (a1 + w*b1) / ((w + 1)*r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function archimedes_arithmetic_mean()
# ----------------------------------------------------------------------
def archimedes_arithmetic_mean(a1, b1, r):
    '''Function Archimedes-Arithmetic-Mean.

         a₁ + b₁
    ac = ───────
           2⋅r
    '''
    # Print method which is used.
    print("*** ARITHMETIC-MEAN ***")
    # Calculate the Archimedes constant.
    ac = (a1 + b1) / (2*r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function archimedes_snellius()
# ----------------------------------------------------------------------
def archimedes_snellius(a1, b1, r):
    '''Function Archimedes-Snellius.

         a₁ + 2⋅b₁
    ac = ─────────
            3⋅r
    '''
    # Print method which is used.
    print("*** SNELLIUS ***")
    # Calculate the Archimedes constant.
    ac = (a1 + 2*b1) / (3*r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function archimedes_doerrie()
# ----------------------------------------------------------------------
def archimedes_doerrie(a1, b1, r):
    '''Function Archimedes-Dörrie.

                        ________
          3⋅a₁⋅b₁    3 ╱      2
         ───────── + ╲╱  a₁⋅b₁
         2⋅a₁ + b₁
    ac = ───────────────────────
               2⋅r
    '''
    # Print method which is used.
    print("*** DÖRRIE ***")
    # Calculate the Archimedes constant.
    ac = ((D(3*a1*b1)/D(2*a1 + b1)) + (D(a1 * b1**2)**(D(1)/D(3))))/D(2*r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function archimedes_netz()
# ----------------------------------------------------------------------
def archimedes_netz(a1, b1, r):
    '''Function Archimedes-Netz.

    Calculate the Archimedes constant using the idea of the so-called
    Snellius acceleration in something like a squared form as well as
    a mixture of harmonic and geometric mean introduced by Dörrie. The
    acceleration of the convergence of the calculation is in principle
    a suitable selected weighted arithmetic mean:
                             _______
               3⋅a₁⋅b₁    3 ╱     2
          4 ⋅ ───────── + ╲╱ a₁⋅b₁
             2⋅a₁ + b₁
     ac = ──────────────────────────
                   5⋅r
    '''
    # Print method which is used.
    print("*** NETZ ***")
    # Calculate the Archimedes constant.
    ac = ((D(12*a1*b1)/D(2*a1 + b1)) + (D(a1 * b1**2)**(D(1)/D(3))))/D(5*r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function archimedes_constant()
# ----------------------------------------------------------------------
def archimedes_constant(a1, b1, r, method=0):
    '''Calculate Pi based on the method.'''
    if method == 0:
        ac = archimedes_netz(a1, b1, r)
    elif method == 1:
        ac = archimedes_doerrie(a1, b1, r)
    elif method == 2:
        ac = archimedes_snellius(a1, b1, r)
    elif method == 3:
        ac = archimedes_arithmetic_mean(a1, b1, r)
    elif method == 4:
        ac = archimedes_weighted_arithmetic_mean(a1, b1, r)
    # Return the Archimedes constant.
    return ac

# ----------------------------------------------------------------------
# Function calculate_pi()
# ----------------------------------------------------------------------
def calculate_pi(iteration=19, r=D(1), method=0, progress=False):
    '''Archimedes algorithm.'''
    # If progress True do something.
    if progress:
        # Hide the cursor.
        hide_cursor()
    # Define the start values.
    a0 = r * 2 * D(3).sqrt()   # half of the outer perimeter
    b0 = r * 3                 # half of the inner perimeter
    # Loop an iteration from 0 to ITERATION plus 1.
    for i in range(0, iteration+1):
        # If progress True do something.
        if progress:
            print_iteration(i)
        # Use the start values in the first loop.
        if i == 0:
            a1 = a0
            b1 = b0
        else:
            # Calculate the half of inner and outer perimeter.
            a1 = D(2*a0*b0)/D(a0 + b0)
            b1 = D(b0*a1).sqrt()
        # Store the old values for the next loop.
        a0 = a1
        b0 = b1
    # Calculate the Archimedes constant using.
    ac = archimedes_constant(a1, b1, r, method=method)
    # If progress True do something.
    if progress:
        # Show the cursor.
        show_cursor()
    # Return the Archimedes constant.
    return str(ac)

# ++++++++++++++++++++
# Main script function
# ++++++++++++++++++++
def main(places, iteration, precision, radius, method, progress, piref):
    '''Main script function.'''
    # Initialise the local variable.
    correct_places = "n/a"
    # Leave script on KeyboardInterrupt exception.
    try:
        # Call the function for calculating Pi.
        ac = calculate_pi(iteration=iteration, r=radius,
                          method=method, progress=progress)
    except KeyboardInterrupt:
        # Clean up and exit script.
        sys.stdout.write("\33[?25h")
        sys.stdout.flush()
        os._exit(1)
    # Print a summary to the screen.
    if PLACES <= len(piref)-2:
        print("\n{0}:".format("Reference"))
        print_pi(str(piref[:places+2]), 50)
        correct_number, correct_places = correct_digits(ac, piref)
    print("{0}:".format("Calculation"))
    print_pi(ac[:places+2], 50)
    print_pi(correct_number[:places+2], 50)
    print("Used precision:", str(precision))
    print("Used iteration:", str(iteration))
    print("\nOutput of requested places:", str(places))
    print("Matching places calculated:", str(correct_places))
    # End of function. Return 1.
    return 1

# Execute script as module or as program.
if __name__ == '__main__':
    # Remove whitespaces from herestring.
    PI = remove_whitespaces(PI_HEREDOC)
    # Call the main script function.
    main(PLACES, ITERATION, PRECISION, RADIUS, METHOD, PROGRESS, PI)
