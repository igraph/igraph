
var svg = d3.select("#igraphx_svg");
var width = svg.attr("width");
var height = svg.attr("height");

svg.append("g").attr("id", "edges");
svg.append("g").attr("id", "vertices");

var color = d3.scale.category20();

var vertices = [];
var edges = [];

var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height])
    .nodes(vertices)
    .links(edges)
    .start();

function findNode(id) {
    for (var i in vertices) {
	if (vertices[i].id == id) { return Number(i); }
    }
    return null;
}

function parseRequest(req) {
    for (var on in req) {	    // on is an operation code, e.g. 'an'
	oper=req[on];		    // oper is the full operation
	for (id in oper) {	    // id is the node or edge id
	    fields=oper[id];        // fields is the hash with the attrs
	    if (on=="an") {
		fields.id=id;
		vertices.push(fields);
	    } else if (on=="ae") {
		fields.id=id;
		fields.source=findNode(fields.source);
		fields.target=findNode(fields.target);
		edges.push(fields);
	    } else if (on=="cn") {
		// TODO
	    } else if (on=="ce") {
		// TODO
	    } else if (on=="dn") {
		// TODO
	    } else if (on=="de") {
		// TODO
	    }
	}
    }
}

function update() {
    var node = svg.select("#vertices").selectAll("circle.node")
	.data(vertices, function(v) { return v.id; });
    
    node.enter().append("circle")
	.attr("class", "node")
	.attr("r", 5)
	.attr("fill", "darkolivegreen")
	.attr("stroke", "darkolivegreen")
	.call(force.drag);
    
    var link = svg.select("#edges").selectAll("line.link")
	.data(edges, function(v) { return v.id; });
        
    link.enter().append("line")
	.attr("class", "link")
	.style("stroke-width", 1)
	.style("stroke", "grey");

    force
	.nodes(vertices)
	.links(edges)
	.start()
	.on("tick", function() {
	    link.attr("x1", function(d) { return d.source.x; })
		.attr("y1", function(d) { return d.source.y; })
		.attr("x2", function(d) { return d.target.x; })
		.attr("y2", function(d) { return d.target.y; });
	    node.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; });
	});
}

var mess;

var socket;
try { 
    socket = new WebSocket("ws://127.0.0.1:7681");

    socket.onopen = function() {
	socket.send("Please!");
    }

    socket.onmessage = function(msg) {
	req=msg.data.split("\n\r");
	for (i in req) {
	    par=JSON.parse(req[i]);
	    parseRequest(par);
	}
	update();
    }
    
    socket.onclose = function() {
    }

} catch(ex) { 
    alert("error: " + ex);
}
