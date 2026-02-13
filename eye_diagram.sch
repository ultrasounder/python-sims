<Qucs Schematic 0.0.19>
<Properties>
  <View=-34,71,571,511,1.85165,0,0>
  <Grid=10,10,1>
  <DataSet=eye_diagram.dat>
  <DataDisplay=eye_diagram.dpl>
  <OpenDisplay=1>
  <Script=eye_diagram.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <SPfile X1 1 180 150 -26 -59 0 0 "C:/Users/vsap7/code/rf-sim/Bad_Via_Model.s2p" 1 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <Vpulse V1 1 20 170 18 -26 0 1 "0 V" 1 "0.9 V" 1 "0" 1 "62.5 ps" 1 "1 ns" 0 "1 ns" 0>
  <Vac V2 1 80 170 18 -26 0 1 "22.5mv" 1 "3210KHz" 0 "0" 0 "0" 0>
  <GND * 1 20 200 0 0 0 0>
  <GND * 1 80 200 0 0 0 0>
  <GND * 1 150 290 0 0 0 0>
  <GND * 1 180 180 0 0 0 0>
  <GND * 1 370 150 0 0 0 0>
  <R R1 1 150 220 15 -26 0 1 "85 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <R R2 1 340 150 -26 15 0 0 "85 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <.TR TR1 1 0 330 0 64 0 0 "lin" 1 "0" 1 "10 us" 1 "100000" 0 "Trapezoidal" 0 "2" 0 "1 ns" 0 "1e-16" 0 "150" 0 "0.001" 0 "1 pA" 0 "1 uV" 0 "26.85" 0 "1e-3" 0 "1e-6" 0 "1" 0 "CroutLU" 0 "no" 0 "yes" 0 "0" 0>
</Components>
<Wires>
  <150 140 150 150 "" 0 0 0 "">
  <80 140 150 140 "" 0 0 0 "">
  <20 140 80 140 "" 0 0 0 "">
  <150 150 150 190 "" 0 0 0 "">
  <150 250 150 290 "" 0 0 0 "">
  <210 150 310 150 "RX_eye" 280 120 40 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
