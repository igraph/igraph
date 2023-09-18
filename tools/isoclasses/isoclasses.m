(* ::Package:: *)

(* ::Text:: *)
(*This Mathematica code is used to generate the "isoclass" lookup tables in the source file isoclasses.c. It is meant to be opened with the Mathematica front end and used as a notebook. However, it is written in plain-text "package" format in order to work well with version control, and to be readable without Mathematica.*)


(* ::Subsubsection:: *)
(*Load IGraph/M*)


Needs["IGraphM`"]


(* ::Subsubsection:: *)
(*Helper functions*)


vec2matD[n_][vec_] := Transpose@MapIndexed[Insert[#1, 0, #2]&, Partition[vec, n-1]]


vec2matUtri[n_][vec_] := PadRight[TakeList[vec, Range[n] - 1], {n, n}]


vec2matU[n_][vec_] :=
  With[{mat=vec2matUtri[n][vec]},
    mat + Transpose[mat]
  ]


(* ::Subsubsection:: *)
(*Directed idx table*)


idxD[n_] := Flatten@vec2matD[n][2^(Range[n (n - 1)] - 1)]


(* ::Subsubsection:: *)
(*Undirected idx table*)


idxU[n_] := Flatten@vec2matU[n][2^(Range[n (n - 1) / 2] - 1)]


(* ::Subsubsection:: *)
(*Directed lookup tables*)


(* ::Text:: *)
(*TODO*)


(* ::Subsubsection:: *)
(*Undirected lookup tables*)


(* ::Text:: *)
(*Vertex count:*)


n = 6;


(* ::Text:: *)
(*How many distinct elements does the adjacency matrix have?*)


k = n (n - 1) / 2


(* ::Subsubsubsection:: *)
(*classedges*)


classedges =
  Module[{mat = vec2matUtri[n][Range[k]], pos},
    pos = Position[mat, _?Positive];
    KeySort@AssociationThread[Extract[mat,pos],pos]
  ] // Values // Flatten;

classedges = Reverse[classedges];


classedges - 1 (* 'classedges' array with 0-based indexing *)


(* ::Subsubsubsection:: *)
(*graphs*)


(* ::Text:: *)
(*Generate the "codes" of all labelled graphs on n vertices.*)


codes = Reverse /@ IntegerDigits[Range[2^k] - 1, 2, k];


graphCodes=
  DeleteDuplicatesBy[
    codes,
    IGBlissCanonicalGraph@*AdjacencyGraph@*vec2matU[n]
  ];


(* 'graphs' array: *)
graphCodesAsInts = FromDigits[#, 2]& /@ Reverse /@ graphCodes


graphs = IGBlissCanonicalGraph @* AdjacencyGraph @* vec2matU[n] /@ graphCodes;


(* 0-based indexes of non-weakly-connected graphs, needed for motif code: *)
Flatten@Position[graphs,g_/;Not@ConnectedGraphQ[g],{1},Heads->False]-1


(* ::Subsubsubsection:: *)
(*lookup table*)


asc = AssociationThread[
  Normal /@ AdjacencyMatrix /@ graphs,
  Range@Length[graphs]
];


lookup = asc @* Normal @* AdjacencyMatrix @* IGBlissCanonicalGraph @* AdjacencyGraph @* vec2matU[n] /@ codes;


(* 'isoclass2' array with 0-based indexing: *)
lookup - 1
