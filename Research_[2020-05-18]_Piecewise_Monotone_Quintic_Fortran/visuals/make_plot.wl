(* Make a plot of the points where the constraints are styled with
   the plot markers (increasing -> up triangle, decreasing -> down
   triangle, extreme point -> circle, flat point -> square) *)
StyledPoints[x_, y_, ms_:100] := {
  ClearAll[dPos, dNeg, dZero, ddZero, left, right,
	red, blue, green, purple, points, series, markers, styles];
  (* Break up all points into increasing, decreasing, extreme, or flat. *)
  dZero={}; ddZero={}; dPos={}; dNeg={};
  For[i=1,i<=Length[x],i++,{
    left=If[i>1,y[[i]]-y[[i-1]],y[[2]]-y[[1]]];
    right=If[i+1<=Length[x],y[[i+1]]-y[[i]],y[[Length[y]]]-y[[Length[y]-1]]];
    Which[
      (Abs[left] < 2^(-26) || Abs[right] < (2^-26)),AppendTo[ddZero,i],
      (left right)<0, AppendTo[dZero,i],
      (left>0),       AppendTo[dPos,i],
      (left<0),       AppendTo[dNeg,i],
      True, Print["ERROR: An error occurred in StyledPoints."]
    ];
  }];
  (* Convert all indices to actual (x,y) points. *)
  toPoints[s=i_]:={x[[i]],y[[i]]};
  (* Set colors for styles. *)
  red = RGBColor[0.84,0.44,0.25];
  blue = RGBColor[0.5,0.76,0.89];
  green = RGBColor[0.26,0.72,0.4];
  purple = RGBColor[0.72,0.47,0.75];
  (* Organize all the points into a list of points, marker styles, and series colors. *)
  series = {}; markers = {}; styles = {};
  If[Length[dZero] > 0, {
    AppendTo[series,  Map[toPoints,dZero]];
    AppendTo[markers, \[FilledSmallCircle]];
    AppendTo[styles,  red];
  }];
  If[Length[ddZero] > 0, {
    AppendTo[series,  Map[toPoints,ddZero]];
    AppendTo[markers, \[FilledSquare]];
    AppendTo[styles,  blue];
  }];
  If[Length[dPos] > 0, {
    AppendTo[series,  Map[toPoints,dPos]];
    AppendTo[markers, \[FilledUpTriangle]];
    AppendTo[styles,  green];
  }];  
  If[Length[dNeg] > 0, {
    AppendTo[series,  Map[toPoints,dNeg]];
    AppendTo[markers, \[FilledDownTriangle]];
    AppendTo[styles,  purple];
  }];  
  (* Create the graphic with all stylized points, return. *)
  points = ListPlot[series[[1]], PlotMarkers->markers, PlotStyle->styles]
}[[1]];


(* Make a plot of the MQSI to data, including stylized markers for all
   data determined by local function change conditions. *)
makePlot[name_, x_, y_, m_:100, width_:400, height_:150] := {
  ClearAll[n, u, qu, data, points, darkGray, graphic];
  Print[""];
  Print["Making a plot for '"<>name<>"'..."];

  (* Read the data that was given to make an MQSI. *)
  n = Length[x];
  u = Subdivide[x[[1]], x[[-1]], m-1];
  
  (* Export a null output file, to make sure one exists. *)
  Export[name<>".out", {0}, "Table", "FieldSeparators" -> " "]; 
  (* Export a data file and points file for the MQSI command line interface. *)
  Export[name<>".data", {n,N[x],N[y]}, "Table", "FieldSeparators" -> " "];
  Export[name<>".pts", {m,N[u]}, "Table", "FieldSeparators" -> " "];

  (* Write a file containing all the points to get the MQSI
     evaluations by using the MQSI command line interface. *)
  Print[""];
  Print["  rm "<>name<>".out"];
  Run["rm "<>name<>".out"]; (* Remove the existing file. *)
  Print["  ./cli "<>name<>".data "<>name<>".pts "<>name<>".out"];
  Run["./cli "<>name<>".data "<>name<>".pts "<>name<>".out"];
  
  (* Read the values of the MQSI at all points. *)
  qu = Map[Internal`StringToDouble,StringSplit[Import[name<>".out","Text"],Whitespace]];

  (* Organize and plot the data. *)
  data = Transpose[{u,qu}];
  Print[""];
  Print["  adding stylized points to graphic.."];
  points = StyledPoints[x,y];
  darkGray = RGBColor[0.34,0.35,0.44];
  Print["  adding MQSI line to graphic.."];
  graphic = Show[
    {ListLinePlot[
      data, PlotStyle->{Thickness[.003],darkGray},
      LabelStyle->{FontSize->10,FontFamily->"Liberation Serif",FontColor->Black}],
     points}, ImageSize->{400,150}, AspectRatio->Full];
  Print["  exporting the graphic to '"<>name<>".eps'.."];
  Export[name<>".eps",graphic];
  Print["  cleaning up files.."];
  Print[""];
  Print["  rm "<>name<>".data "<>name<>".pts "<>name<>".out "];
  Run["rm "<>name<>".data "<>name<>".pts "<>name<>".out "]; (* Remove the files *)
}

makePlot["testing", {1, 2, 3, 4}, {0, 1, 1.1, 2}]

