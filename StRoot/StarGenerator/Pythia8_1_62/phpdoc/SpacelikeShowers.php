<html>
<head>
<title>Spacelike Showers</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='SpacelikeShowers.php'>

<h2>Spacelike Showers</h2>

The PYTHIA algorithm for spacelike initial-state showers is 
based on the article [<a href="Bibliography.php" target="page">Sjo05</a>], where a 
transverse-momentum-ordered backwards evolution scheme is introduced, 
with the extension to fully interleaved evolution covered in 
[<a href="Bibliography.php" target="page">Cor10a</a>]. 
This algorithm is a further development of the virtuality-ordered one 
presented in [<a href="Bibliography.php" target="page">Sj085</a>], with matching to first-order matrix 
element for <i>Z^0</i>, <i>W^+-</i> and Higgs (in the 
<i>m_t -> infinity</i> limit) production as introduced in 
[<a href="Bibliography.php" target="page">Miu99</a>]. 

<p/>
The normal user is not expected to call <code>SpaceShower</code> 
directly, but only have it called from <code>Pythia</code>, 
via <code>PartonLevel</code>. Some of the parameters below, 
in particular <code>SpaceShower:alphaSvalue</code>, 
would be of interest for a tuning exercise, however. 

<h3>Main variables</h3>

The maximum <i>pT</i> to be allowed in the shower evolution is
related to the nature of the hard process itself. It involves a
delicate balance between not doublecounting and not leaving any
gaps in the coverage. The best procedure may depend on information 
only the user has: how the events were generated and mixed (e.g. with 
Les Houches Accord external input), and how they are intended to be 
used. Therefore a few options are available, with a sensible default 
behaviour.

<br/><br/><table><tr><td><strong>SpaceShower:pTmaxMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Way in which the maximum shower evolution scale is set to match the 
scale of the hard process itself.
<br/>
<input type="radio" name="1" value="0" checked="checked"><strong>0 </strong>: <b>(i)</b> if the final state of the hard process  (not counting subsequent resonance decays) contains at least one quark  (<ei>u, d, s, c ,b</ei>), gluon or photon then <ei>pT_max</ei>  is chosen to be the factorization scale for internal processes  and the <code>scale</code> value for Les Houches input;  <b>(ii)</b> if not, emissions are allowed to go all the way up to  the kinematical limit.  The reasoning is that in the former set of processes the ISR emission of yet another quark, gluon or photon could lead to doublecounting, while no such danger exists in the latter case. <br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: always use the factorization scale for an internal process and the <code>scale</code> value for Les Houches input,  i.e. the lower value. This should avoid doublecounting, but may leave out some emissions that ought to have been simulated. (Also known as wimpy showers.) <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: always allow emissions up to the kinematical limit. This will simulate all possible event topologies, but may lead to doublecounting.  (Also known as power showers.) <br/>
<br/><b>Note 1:</b> These options only apply to the hard interaction.
Emissions off subsequent multiparton interactions are always constrainted
to be below the factorization scale of the process itself.  
<br/><b>Note 2:</b> Some processes contain matrix-element matching
to the first emission; this is the case notably for single 
<ei>gamma^*/Z^0, W^+-</ei> and <ei>H^0</ei> production. Then default
and option 2 give the correct result, while option 1 should never
be used. 

<br/><br/><table><tr><td><strong>SpaceShower:pTmaxFudge </td><td></td><td> <input type="text" name="2" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
In cases where the above <code>pTmaxMatch</code> rules would imply
that <i>pT_max = pT_factorization</i>, <code>pTmaxFudge</code> 
introduces a multiplicative factor <i>f</i> such that instead
<i>pT_max = f * pT_factorization</i>. Only applies to the hardest
interaction in an event, cf. below. It is strongly suggested that 
<i>f = 1</i>, but variations around this default can be useful to 
test this assumption. 
  

<br/><br/><table><tr><td><strong>SpaceShower:pTmaxFudgeMPI </td><td></td><td> <input type="text" name="3" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 2.0</code>)</td></tr></table>
A multiplicative factor <i>f</i> such that 
<i>pT_max = f * pT_factorization</i>, as above, but here for the
non-hardest interactions (when multiparton interactions are allowed).
  

<br/><br/><table><tr><td><strong>SpaceShower:pTdampMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
These options only take effect when a process is allowed to radiate up 
to the kinematical limit by the above <code>pTmaxMatch</code> choice, 
and no matrix-element corrections are available. Then, in many processes,
the fall-off in <ei>pT</ei> will be too slow by one factor of <ei>pT^2</ei>. 
That is, while showers have an approximate <ei>dpT^2/pT^2</ei> shape, often 
it should become more like <ei>dpT^2/pT^4</ei> at <ei>pT</ei> values above
the scale of the hard process. Whether this actually is the case 
depends on the particular process studied, e.g. if <ei>t</ei>-channel 
gluon exchange is likely to dominate. If so, the options below could
provide a reasonable high-<ei>pT</ei> behaviour without requiring 
higher-order calculations. 
<br/>
<input type="radio" name="4" value="0" checked="checked"><strong>0 </strong>: emissions go up to the kinematical limit,  with no special dampening. <br/>
<input type="radio" name="4" value="1"><strong>1 </strong>: emissions go up to the kinematical limit,   but dampened by a factor <ei>k^2 Q^2_fac/(pT^2 + k^2 Q^2_fac)</ei>,  where <ei>Q_fac</ei> is the factorization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below. <br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: emissions go up to the kinematical limit,  but dampened by a factor <ei>k^2 Q^2_ren/(pT^2 + k^2 Q^2_ren)</ei>,  where <ei>Q_ren</ei> is the renormalization scale and <ei>k</ei> is a  multiplicative fudge factor stored in <code>pTdampFudge</code> below.  <br/>
<br/><b>Note:</b> These options only apply to the hard interaction.
Emissions off subsequent multiparton interactions are always constrainted
to be below the factorization scale of the process itself.  

<br/><br/><table><tr><td><strong>SpaceShower:pTdampFudge </td><td></td><td> <input type="text" name="5" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.25</code>; <code>maximum = 4.0</code>)</td></tr></table>
In cases 1 and 2 above, where a dampening is imposed at around the
factorization or renormalization scale, respectively, this allows the
<i>pT</i> scale of dampening of radiation by a half to be shifted 
by this factor relative to the default <i>Q_fac</i> or <i>Q_ren</i>. 
This number ought to be in the neighbourhood of unity, but variations 
away from this value could do better in some processes.
  

<p/>
The amount of QCD radiation in the shower is determined by 
<br/><br/><table><tr><td><strong>SpaceShower:alphaSvalue </td><td></td><td> <input type="text" name="6" value="0.137" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.137</strong></code>; <code>minimum = 0.06</code>; <code>maximum = 0.25</code>)</td></tr></table>
The <i>alpha_strong</i> value at scale <code>M_Z^2</code>. 
Default value is picked equal to the one used in CTEQ 5L.  
  

<p/>
The actual value is then regulated by the running to the scale 
<i>pT^2</i>, at which it is evaluated
<br/><br/><table><tr><td><strong>SpaceShower:alphaSorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Order at which <ei>alpha_strong</ei> runs,
<br/>
<input type="radio" name="7" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_strong</ei> is kept  fixed.<br/>
<input type="radio" name="7" value="1" checked="checked"><strong>1 </strong>: first order, which is the normal value.<br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: second order. Since other parts of the code do  not go to second order there is no strong reason to use this option,  but there is also nothing wrong with it.<br/>

<p/>
QED radiation is regulated by the <i>alpha_electromagnetic</i>
value at the <i>pT^2</i> scale of a branching.
 
<br/><br/><table><tr><td><strong>SpaceShower:alphaEMorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
The running of <ei>alpha_em</ei>.
<br/>
<input type="radio" name="8" value="1" checked="checked"><strong>1 </strong>: first-order running, constrained to agree with <code>StandardModel:alphaEMmZ</code> at the <ei>Z^0</ei> mass. <br/>
<input type="radio" name="8" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed at its value at vanishing momentum transfer.<br/>
<input type="radio" name="8" value="-1"><strong>-1 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept  fixed, but at <code>StandardModel:alphaEMmZ</code>, i.e. its value at the <ei>Z^0</ei> mass. <br/>

<p/>
There are two complementary ways of regularizing the small-<i>pT</i> 
divergence, a sharp cutoff and a smooth dampening. These can be 
combined as desired but it makes sense to coordinate with how the 
same issue is handled in multiparton interactions.

<br/><br/><strong>SpaceShower:samePTasMPI</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Regularize the <i>pT -> 0</i> divergence using the same sharp cutoff 
and smooth dampening parameters as used to describe multiparton interactions.
That is, the <code>MultipartonInteractions:pT0Ref</code>, 
<code>MultipartonInteractions:ecmRef</code>, 
<code>MultipartonInteractions:ecmPow</code> and 
<code>MultipartonInteractions:pTmin</code> parameters are used to regularize 
all ISR QCD radiation, rather than the corresponding parameters below. 
This is a sensible physics ansatz, based on the assumption that colour 
screening effects influence both MPI and ISR in the same way. Photon 
radiation is regularized separately in either case.
<br/><b>Warning:</b> if a large <code>pT0</code> is picked for multiparton 
interactions, such that the integrated interaction cross section is 
below the nondiffractive inelastic one, this <code>pT0</code> will 
automatically be scaled down to cope. Information on such a rescaling 
does NOT propagate to <code>SpaceShower</code>, however.
   

<p/>
The actual <code>pT0</code> parameter used at a given CM energy scale, 
<i>ecmNow</i>, is obtained as
<br/><i>
    pT0 = pT0(ecmNow) = pT0Ref * (ecmNow / ecmRef)^ecmPow 
</i><br/>
where <i>pT0Ref</i>, <i>ecmRef</i> and <i>ecmPow</i> are the 
three parameters below.

<br/><br/><table><tr><td><strong>SpaceShower:pT0Ref </td><td></td><td> <input type="text" name="10" value="2.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.0</strong></code>; <code>minimum = 0.5</code>; <code>maximum = 10.0</code>)</td></tr></table>
Regularization of the divergence of the QCD emission probability for 
<i>pT -> 0</i> is obtained by a factor <i>pT^2 / (pT0^2 + pT^2)</i>, 
and by using an <i>alpha_s(pT0^2 + pT^2)</i>. An energy dependence 
of the <i>pT0</i> choice is introduced by the next two parameters, 
so that <i>pT0Ref</i> is the <i>pT0</i> value for the reference 
cm energy, <i>pT0Ref = pT0(ecmRef)</i>.   
  

<br/><br/><table><tr><td><strong>SpaceShower:ecmRef </td><td></td><td> <input type="text" name="11" value="1800.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1800.0</strong></code>; <code>minimum = 1.</code>)</td></tr></table>
The <i>ecmRef</i> reference energy scale introduced above.
  

<br/><br/><table><tr><td><strong>SpaceShower:ecmPow </td><td></td><td> <input type="text" name="12" value="0.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0</strong></code>; <code>minimum = 0.</code>; <code>maximum = 0.5</code>)</td></tr></table>
The <i>ecmPow</i> energy rescaling pace introduced above.
  

<br/><br/><table><tr><td><strong>SpaceShower:pTmin </td><td></td><td> <input type="text" name="13" value="0.2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.2</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.0</code>)</td></tr></table>
Lower cutoff in <i>pT</i>, below which no further ISR branchings 
are allowed. Normally the <i>pT0</i> above would be used to 
provide the main regularization of the branching rate for 
<i>pT -> 0</i>, in which case <i>pTmin</i> is used  mainly for 
technical reasons. It is possible, however, to set <i>pT0Ref = 0</i> 
and use <i>pTmin</i> to provide a step-function regularization, 
or to combine them in intermediate approaches. Currently <i>pTmin</i> 
is taken to be energy-independent.  
  

<br/><br/><table><tr><td><strong>SpaceShower:pTminChgQ </td><td></td><td> <input type="text" name="14" value="0.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.5</strong></code>; <code>minimum = 0.01</code>)</td></tr></table>
Parton shower cut-off <i>pT</i> for photon coupling to a coloured 
particle.
  

<br/><br/><table><tr><td><strong>SpaceShower:pTminChgL </td><td></td><td> <input type="text" name="15" value="0.0005" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.0005</strong></code>; <code>minimum = 0.0001</code>)</td></tr></table>
Parton shower cut-off mass for pure QED branchings. 
Assumed smaller than (or equal to) <i>pTminChgQ</i>.
  

<br/><br/><strong>SpaceShower:rapidityOrder</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Force emissions, after the first,  to be ordered in rapidity,
i.e. in terms of decreasing angles in a backwards-evolution sense. 
Could be used to probe sensitivity to unordered emissions.
Only affects QCD emissions.
  

<h3>Further variables</h3>

These should normally not be touched. Their only function is for
cross-checks.

<p/>
There are three flags you can use to switch on or off selected
branchings in the shower: 

<br/><br/><strong>SpaceShower:QCDshower</strong>  <input type="radio" name="17" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="17" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow a QCD shower; on/off = true/false.
  

<br/><br/><strong>SpaceShower:QEDshowerByQ</strong>  <input type="radio" name="18" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="18" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow quarks to radiate photons; on/off = true/false.
  

<br/><br/><strong>SpaceShower:QEDshowerByL</strong>  <input type="radio" name="19" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="19" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Allow leptons to radiate photons; on/off = true/false.
  

<p/>
There are some further possibilities to modify the shower:

<br/><br/><strong>SpaceShower:MEcorrections</strong>  <input type="radio" name="20" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="20" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections; on/off = true/false.
  

<br/><br/><strong>SpaceShower:MEafterFirst</strong>  <input type="radio" name="21" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="21" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use of matrix element corrections also after the first emission,
for dipole ends of the same system that did not yet radiate.
Only has a meaning if <code>MEcorrections</code> above is 
switched on. 
  

<br/><br/><strong>SpaceShower:phiPolAsym</strong>  <input type="radio" name="22" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="22" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Azimuthal asymmetry induced by gluon polarization; on/off = true/false.
  

<br/><br/><strong>SpaceShower:phiIntAsym</strong>  <input type="radio" name="23" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="23" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Azimuthal asymmetry induced by interference; on/off = true/false.
  

<br/><br/><table><tr><td><strong>SpaceShower:strengthIntAsym </td><td></td><td> <input type="text" name="24" value="0.7" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.7</strong></code>; <code>minimum = 0.</code>; <code>maximum = 0.9</code>)</td></tr></table>
Size of asymmetry induced by interference. Natural value of order 0.5; 
expression would blow up for a value of 1.
  

<br/><br/><table><tr><td><strong>SpaceShower:nQuarkIn  </td><td></td><td> <input type="text" name="25" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed quark flavours in <i>g -> q qbar</i> branchings,
when kinematically allowed, and thereby also in incoming beams. 
Changing it to 4 would forbid <i>g -> b bbar</i>, etc.
  

<h3>Technical notes</h3>

Almost everything is equivalent to the algorithm in [1]. Minor changes 
are as follows.
<ul>
<li>
It is now possible to have a second-order running <i>alpha_s</i>,
in addition to fixed or first-order running. 
</li>
<li>
The description of heavy flavour production in the threshold region 
has been modified, so as to be more forgiving about mismatches 
between the <i>c/b</i>  masses used in Pythia relative to those 
used in a respective PDF parametrization. The basic idea is that, 
in the threshold region of a heavy quark <i>Q</i>, <i>Q = c/b</i>, 
the effect of subsequent <i>Q -> Q g</i> branchings is negligible. 
If so, then
<br/><i>
   f_Q(x, pT2) = integral_mQ2^pT2  dpT'2/pT'2 * alpha_s(pT'2)/2pi
      * integral P(z) g(x', pT'2) delta(x - z x')
</i><br/>
so use this to select the <i>pT2</i> of the <i>g -> Q Qbar</i> 
branching. In the old formalism the same kind of behaviour should 
be obtained, but by a cancellation of a <i>1/f_Q</i> that diverges 
at the theshold and a Sudakov that vanishes.
<br/>
The strategy therefore is that, once <i>pT2 &lt; f * mQ2</i>, with 
<i>f</i> a parameter of the order of 2, a <i>pT2</i> is chosen 
like <i>dpT2/pT2</i> between <i>mQ2</i> and <i>f * mQ2</i>, a
nd a <i>z</i> flat in the allowed range. Thereafter acceptance
is based on the product of three factors, representing the running
of <i>alpha_strong</i>, the splitting kernel (including the mass term) 
and the gluon density weight. At failure, a new <i>pT2</i> is chosen 
in the same  range, i.e. is not required to be lower since no Sudakov 
is involved. 
</li>
<li>
The QED algorithm now allows for hadron beams with non-zero photon
content. The backwards-evolution of a photon in a hadron is identical
to that of a gluon, with <i>CF -> eq^2</i> and <i>CA -> 0</i>.
Note that this will only work in conjunction with 
parton distribution that explicitly include photons as part of the
hadron structure (such as the MRST2004qed set). Since Pythia's
internal sets do not allow for photon content in hadrons, it is thus 
necessary to use the LHAPDF interface to make use of this feature. The
possibility of a fermion backwards-evolving to a photon has not yet
been included, nor has photon backwards-evolution in lepton beams. 
</li>
</ul>

<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "0")
{
$data = "SpaceShower:pTmaxMatch = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "1.0")
{
$data = "SpaceShower:pTmaxFudge = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1.0")
{
$data = "SpaceShower:pTmaxFudgeMPI = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "SpaceShower:pTdampMatch = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1.0")
{
$data = "SpaceShower:pTdampFudge = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0.137")
{
$data = "SpaceShower:alphaSvalue = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1")
{
$data = "SpaceShower:alphaSorder = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1")
{
$data = "SpaceShower:alphaEMorder = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "SpaceShower:samePTasMPI = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "2.0")
{
$data = "SpaceShower:pT0Ref = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1800.0")
{
$data = "SpaceShower:ecmRef = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0.0")
{
$data = "SpaceShower:ecmPow = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "0.2")
{
$data = "SpaceShower:pTmin = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0.5")
{
$data = "SpaceShower:pTminChgQ = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "0.0005")
{
$data = "SpaceShower:pTminChgL = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "SpaceShower:rapidityOrder = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "on")
{
$data = "SpaceShower:QCDshower = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "on")
{
$data = "SpaceShower:QEDshowerByQ = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "on")
{
$data = "SpaceShower:QEDshowerByL = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "on")
{
$data = "SpaceShower:MEcorrections = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "on")
{
$data = "SpaceShower:MEafterFirst = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "on")
{
$data = "SpaceShower:phiPolAsym = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "on")
{
$data = "SpaceShower:phiIntAsym = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0.7")
{
$data = "SpaceShower:strengthIntAsym = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "5")
{
$data = "SpaceShower:nQuarkIn = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2012 Torbjorn Sjostrand -->

