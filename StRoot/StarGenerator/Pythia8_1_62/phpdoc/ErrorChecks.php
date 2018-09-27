<html>
<head>
<title>Error Checks</title>
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

<form method='post' action='ErrorChecks.php'>

<h2>Error Checks</h2>

There is a few settings related to error checking during program
execution. Many other checks are performed as well, but do not 
have any specific setting related to themselves.

<br/><br/><strong>Check:abortIfVeto</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
There are a few ways in which an event can be vetoed, the most
common being a <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>User Hooks</a> test.
Normally this will simply mean that the next parton-level
configuration is selected inside the <code>Pythia::next()</code>
routine, without any need for a user intervention. With this
option switched on, however, <code>Pythia::next()</code> will 
return <code>false</code>. It is then up to the user to decide
what to do next. 
  

<br/><br/><strong>Check:particleData</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Check the particle data tables for potential problems during 
initialization. This includes inconsistent use of charge in particle 
names, inconsistent setup of mass, mass range, width and lifetime, 
sum of branching ratios not unity (allowed but discouraged) or charge 
not conserved in a decay channel. Warnings should be viewed as reasons 
to check further, but need not indicate a true problem, and also not all 
problems may be caught. 
The <code>pythia.particleData.checkTable(level)</code> method,
used for these checks, may also be called directly. 
  

<br/><br/><table><tr><td><strong>Check:levelParticleData  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>)</td></tr></table>
The level of verbosity and checks of particle data, if switched on.
<br/>
<input type="radio" name="3" value="0"><strong>0 </strong>: mimimal amount of checks, e.g. that no channels open. <br/>
<input type="radio" name="3" value="1" checked="checked"><strong>1 </strong>: further warning if individual channels closed, except for resonances.<br/>
<input type="radio" name="3" value="2"><strong>2 </strong>: also print branching-ratio-averaged threshold mass except for resonances.<br/>
<input type="radio" name="3" value="11"><strong>11 </strong>: as 1, but include resonances in detailed checks. <br/>
<input type="radio" name="3" value="12"><strong>12 </strong>: as 2, but include resonances in detailed checks. <br/>

<br/><br/><strong>Check:event</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When an event has been successfully generated, check that the 
final event record in <code>event</code> does not contain any 
unphysical particles, or nonconserved charge or energy-momentum. 
If this check fails, then <code>pythia.next()</code> obtains the 
value <code>false</code>, i.e. the event is counted as aborted.
  

<br/><br/><table><tr><td><strong>Check:nErrList  </td><td></td><td> <input type="text" name="5" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
The number of erroneous events, in the above check, for which 
event listing and other detailed information will be printed. 
After that, only the normal error messages will be issued. 
Error counters are always updated, and accumulated numbers can be   
shown with <code>pythia.statistics()</code> at the end of the run.
  

<br/><br/><table><tr><td><strong>Check:epTolErr </td><td></td><td> <input type="text" name="6" value="1e-4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-4</strong></code>)</td></tr></table>
Maximum allowed summed deviation of <i>E</i>, <i>p_x</i>, 
<i>p_y</i> and <i>p_z</i> between the incoming beams and the 
final state, as a fraction of the initial energy, above which the 
event is counted as aborted.
(Unfortunetely roundoff errors do not scale linearly with the energy, 
and also have a very long tail. So while most events at lower energies 
may be correct to better than 1e-10, at LHC it does not have to signal 
any fundamental bug if also the default tolerance above is violated 
occasionally.)
  

<br/><br/><table><tr><td><strong>Check:epTolWarn </td><td></td><td> <input type="text" name="7" value="1e-6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-6</strong></code>)</td></tr></table>
A check on the same summed deviation as above, but counted as a 
warning rather than an error, and not leading to the event being
classified as aborted.
  

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

if($_POST["1"] != "off")
{
$data = "Check:abortIfVeto = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Check:particleData = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "1")
{
$data = "Check:levelParticleData = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "Check:event = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "0")
{
$data = "Check:nErrList = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1e-4")
{
$data = "Check:epTolErr = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1e-6")
{
$data = "Check:epTolWarn = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2012 Torbjorn Sjostrand -->
