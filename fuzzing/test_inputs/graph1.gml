Creator	"yFiles"
Version	"2.18"
graph
[
	hierarchic	1
	label	""
	directed	1
	node
	[
		id	0
		label	"A"
		graphics
		[
			x	15.0
			y	125.48405928593465
			w	30.0
			h	30.0
			type	"ellipse"
			raisedBorder	0
			fill	"#FFCC00"
			outline	"#000000"
		]
		LabelGraphics
		[
			text	"A"
			fontSize	12
			fontName	"Dialog"
			model	"null"
		]
	]
	node
	[
		id	1
		label	"B"
		graphics
		[
			x	59.830739982104404
			y	70.46602381781756
			w	30.0
			h	30.0
			type	"octagon"
			raisedBorder	0
			fill	"#33CCCC"
			outline	"#000000"
		]
		LabelGraphics
		[
			text	"B"
			fontSize	12
			fontName	"Dialog"
			model	"null"
		]
	]
	node
	[
		id	2
		label	"C"
		graphics
		[
			x	104.94041508180403
			y	15.0
			w	30.0
			h	30.0
			type	"ellipse"
			raisedBorder	0
			fill	"#FF9900"
			outline	"#000000"
		]
		LabelGraphics
		[
			text	"C"
			fontSize	12
			fontName	"Dialog"
			model	"null"
		]
	]
	edge
	[
		source	0
		target	1
		label	"edge 2"
		value 42
		graphics
		[
			type	"arc"
			fill	"#000000"
			targetArrow	"standard"
			arcType	"fixedRatio"
			arcHeight	17.74256706237793
			arcRatio	1.0
			Line
			[
				point
				[
					x	15.0
					y	125.48405928593465
				]
				point
				[
					x	19.14604949951172
					y	86.7673568725586
				]
				point
				[
					x	59.830739982104404
					y	70.46602381781756
				]
			]
		]
		edgeAnchor
		[
			xSource	-0.3009873023843175
			xTarget	-0.3009873023843183
		]
		LabelGraphics
		[
			text	"edge 2"
			fontSize	12
			fontName	"Dialog"
			configuration	"AutoFlippingLabel"
			contentWidth	43.791015625
			contentHeight	18.1328125
			model	"null"
			position	"null"
		]
	]
	edge
	[
		source	1
		target	0
		label	"edge 3"
		graphics
		[
			type	"arc"
			fill	"#000000"
			targetArrow	"standard"
			arcType	"fixedRatio"
			arcHeight	17.742568969726562
			arcRatio	1.0
			Line
			[
				point
				[
					x	59.830739982104404
					y	70.46602381781756
				]
				point
				[
					x	55.684688568115234
					y	109.18273162841797
				]
				point
				[
					x	15.0
					y	125.48405928593465
				]
			]
		]
		edgeAnchor
		[
			xSource	0.30098730238431803
			xTarget	0.3009873023843179
		]
		value 1.23
		LabelGraphics
		[
			text	"edge 3"
			fontSize	12
			fontName	"Dialog"
			configuration	"AutoFlippingLabel"
			contentWidth	43.791015625
			contentHeight	18.1328125
			model	"null"
			position	"null"
		]
	]
	edge
	[
		source	1
		target	2
		label	"edge 1"
		graphics
		[
			type	"arc"
			fill	"#FF0000"
			targetArrow	"standard"
			arcType	"fixedRatio"
			arcHeight	17.87344741821289
			arcRatio	1.0
			Line
			[
				point
				[
					x	59.830739982104404
					y	70.46602381781756
				]
				point
				[
					x	68.5190658569336
					y	31.455595016479492
				]
				point
				[
					x	104.94041508180403
					y	15.0
				]
			]
		]
		LabelGraphics
		[
			text	"edge 1"
			fontSize	12
			fontName	"Dialog"
			configuration	"AutoFlippingLabel"
			contentWidth	43.791015625
			contentHeight	18.1328125
			model	"null"
			position	"null"
		]
		value 34.27
	]
]
