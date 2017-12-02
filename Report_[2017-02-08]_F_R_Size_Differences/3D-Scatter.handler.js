var graph3d_array = graph3d_array || [];

function redraw_plots() {
    for (var i = 0; i < graph3d_array.length; i++) {
	graph3d_array[i].redraw();
    }
}

// Redraw the plots if the window changes size
window.addEventListener('onresize', redraw_plots);

function add_3D_Scatter(group, project, view, name, selected_data) {
// Function for creating 3D plot based on JSON returned by API call.
    function create_options (x_axis, y_axis, z_axis) {
    // Generate an "options" object for visualizing the 3D scatter plot
	var options = {	
	    align: 'center',
	    width:  '100%',
	    height: '60vmin',
	    style: "dot",
	    backgroundColor: {fill:'transparent', stroke:'transparent', strokWidth:0},
	    cameraPosition: {horizontal: 1.0, vertical: 0.5, distance: 2.3},
            dotSizeRatio: 0.015,
            keepAspectRatio: false,
            // verticalAspectRatio: 0.8,
	    // showShadow: false,
            xValueLabel: function (value) {
		return (value.toExponential(1));
            },
            yValueLabel: function (value) {
		return (value.toExponential(1));
            },
            zValueLabel: function (value) {
		return (value.toExponential(1));
            },
            xLabel: x_axis,
            yLabel: y_axis,
            zLabel: z_axis,
	    // Option tooltip can be true, false, or a function returning a string with HTML contents
            tooltip: function (point) {
		// parameter point contains properties x, y, z
		return ('('+(point.x).toFixed(0)+','+
			(point.y).toFixed(0)+'):  '+
			(point.z).toExponential(1)) ;
            }
	};
	return options;
    };

    // Beginning of 'add_3D_Scatter' function
    var plot_id = name + "--" + "3D-Scatter"
    var plot_title = document.createElement("div");
    var plot = document.createElement("div");
    var holder = document.createElement("div");

    plot_title.id = plot_id + "--title";
    plot.id = plot_id + "--plot";
    plot_title.style.paddingTop = "10px";
    holder.style.backgroundColor = "#fafafa";
    // holder.style.borderWidth =  "1px";
    // holder.style.borderColor = "#e5e5e5";
    // holder.style.borderStyle = "solid";
    // holder.style.borderRadius = "10px";
    holder.onshow = function () { redraw_plots() };
    holder.onload = function () {console.log("ERROR: 'onload' called before ready.")};
    holder.appendChild(plot_title);
    holder.appendChild(plot);

    // Make the request that will get the view for this 3D figure
    var request = $.getJSON("/View/"+group+"/"+project+"/3D-Scatter/"+
			    [view,name,selected_data].join("/"));
    
    // Handle JSON returned by server
    request.then(function (resp) {
	// resp = { x_axis:<string>,
	//          y_axis:<string>,
	//          z_axis:<string>,
	//          points:<array of <array[3] of <float>>>
	// }

	// Set options for viewing the 3D scatter plot
	var options = create_options(resp.x_axis, resp.y_axis, resp.z_axis);
	// Create and populate a data set for this visualization.
	var data = new vis.DataSet();
	for (var i = 0; i < resp.points.length; i+=1) {
	    data.add({x:resp.points[i][0], y:resp.points[i][1], z:resp.points[i][2]})
	}   
	// Create the visualization
	var graph3d = new vis.Graph3d(plot, data, options);
	graph3d_array.push(graph3d);

	// Insert a title bar before the 3D visualization and set its contents
	plot_title.innerHTML = replace_all(
	    "<b>"+view+":</b> " + name, replace.old, replace.new);

	// Call 'done' function that should have been overwritten by parent
	holder.onload()	
    });

    // Return the holder that has now been almost fully initialized,
    // it is expected that holder.onload will be overwritten.
    return holder;
}

var handlers = handlers || {};
handlers["3D-Scatter"] = add_3D_Scatter;
