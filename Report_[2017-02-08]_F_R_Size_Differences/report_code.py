# Global variables
import sys
import os
home = os.path.expanduser("~")
sys.path = [home + "/Git_Analytics/Analytics-Vis/Code"] + sys.path
curr_dir = os.path.abspath(os.path.dirname(__file__))

from settings import data_for_spec, spec_str_to_dict
import pygal
import numpy as np

NUM_POINTS = 500
CUSTOM_CSS_FILE = 'file://'+home+'/Git_Analytics/Analytics-Vis/Code/ViewData/Value_Distribution/no_dot.css'

# Takes: .format(title, mean_diff, stdev_diff, var_diff, fig_name)
REPORT_ENTRY = '''
    <div class="entry">
      <div class="entry-title" onclick="$(this.parentNode.children[1]).collapse('toggle')">
	%s
      </div>
      <div class="holder collapse in" aria-expanded="true">
	<div class="entry-text">
	  <p style="text-decoration: underline; font-size: 16px; color: #333;">Relative Increase from Min to Max</p>
	  <strong>Mean:</strong> %.2f%%
	  <br>
	  <strong>Stdev:</strong> %.2f%%
	  <br>
	  <strong>Variance:</strong> %.2f%%
	</div>
	<figure>
	  <embed type="image/svg+xml" src="./%s"> </embed>
	</figure>
        <br style="clear: both">
      </div>      
    </div>
'''


VAR_PLOT = '''
      var points = %s;
      // Set options for viewing the 3D scatter plot
      var options = create_options("%s", "%s", "%s");
      // Create and populate a data set for this visualization.
      var data = new vis.DataSet();
      for (var i = 0; i < points.length; i+=1) {
          data.add({x:points[i][0], y:points[i][1], z:points[i][2]})
      }   
      // Create the visualization
      var graph3d = new vis.Graph3d(plot, data, options);
'''

HTML_BEGIN = '''
<!DOCTYPE html>
<html lang="en">
  <head>
    <title>Report on F Size and R Size Variance</title>
    <!-- Default setup for a webpage with basic bootstrap -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="./bootstrap.min.css">
    <link rel="stylesheet" href="./analysis.css">

    <script src="./vis.js"></script>
    <script src="./jquery.min.js"></script>
    <script src="./bootstrap.min.js"></script>
    <script src="./pygal-tooltips.min.js"></script>
    <script src="./3D-Scatter.handler.js"></script>

    <style type="text/css">
      #Brand { cursor: default; }
      #Project { float: right; }

      #Project-Name:hover { color: #555; }
      #Project-Name {
      float: right;
      font-size: 16px;
      padding: 14px 5px;
      color: #777;
      text-decoration: none;
      cursor: default;
      }

      div.entry {      
      border-top: 1px solid #eee;
      border-bottom: 1px solid #eee;
      background-color: #fdfdfd;
      width: 90vw;
      margin: auto;
      margin-top: 10px;
      margin-bottom: 10px;
      }
      div.entry-title {
      padding: 5px 0px;
      font-size: 14px;
      }
      div.holder {
      padding: 10px;
      }
      div.entry-text {
      font-size: 14px;
      float: left;
      display: inline-block;
      padding-top: 15vw;
      width: 25vw;
      height: 45vw;
      }
      figure {
      width: 60vw;
      float: right;
      display: inline-block;
      }

    </style>
  </head>

  <body id="Body">
    <!-- Nav bar -->
    <div class="nav">
	<a id="Brand">Virginia Tech</a>
	<div id="Project">
	  <a id="Project-Name">VarSys</a>
	</div>
    </div>

    <div class="entry">
      <div class="entry-title" onclick="$(this.parentNode.children[1]).collapse('toggle')">
	<strong>Scatter Plot of Variance for 2 Independent 40 Runs</strong>
      </div>
      <div id="Variance-Scatter" class="collapse in" aria-expanded="true">
      </div>
    </div>

'''
HTML_MID = '''
    <script>
      function create_options (x_axis, y_axis, z_axis) {
      // Generate an "options" object for visualizing the 3D scatter plot
    	var options = {
    	    align: 'center',
    	    width:  '90%',
    	    height: '60vmin',
    	    style: "dot",
    	    backgroundColor: {fill:'transparent', stroke:'transparent', strokWidth:0},
    	    cameraPosition: {horizontal: 1.0, vertical: 0.5, distance: 2.8},
            dotSizeRatio: 0.015,
            keepAspectRatio: false,
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

      var plot = document.getElementById("Variance-Scatter");
'''
HTML_END = '''
    </script>
  </body>
</html>
'''

# # Make a plot that stacks, only displaying where data was overlapping and non-overlapping
# def disparity_plot(data1, data2, file_name):
#     data1.sort()
#     data2.sort()
#     min_pt = min(min(data1), min(data2))
#     max_pt = max(max(data1), max(data2))

#     agreement = []
#     disparity = []
#     for i in range(NUM_POINTS):
#         value = (max_pt - min_pt) * (i / NUM_POINTS) + min_pt
#         v1 = sum(data1 < value) / len(data1)
#         v2 = sum(data2 < value) / len(data2)
#         agreement.append(min(v1,v2))
#         disparity.append(max(v1,v2) - min(v1,v2))

#     # Create the plot
#     config = pygal.Config(fill=True, show_legend=True)
#     print("Loading custom CSS file:",CUSTOM_CSS_FILE)
#     config.css.append(CUSTOM_CSS_FILE)
#     chart = pygal.StackedLine(config)

#     chart.title = "Cumulative Distribution function"
#     chart.x_title = "Throughput"
#     chart.y_title = "Percentage of data"
#     chart.x_label_rotation = 10
#     chart.add("Agreement :)", agreement)
#     chart.add("Disparity :(", disparity)
#     chart.render_to_file(file_name)    

# Make a CDF plot (as a line)
def cdf_plot(data1, data2, file_name):
    data1.sort()
    data2.sort()
    min_pt = min(min(data1), min(data2))
    max_pt = max(max(data1), max(data2))

    data1_pts = [(min_pt,0)]
    data2_pts = [(min_pt,0)]
    for i,(val1, val2) in enumerate(zip(data1, data2)):
        data1_pts.append((val1, 100.0*i/(len(data1)-1)))
        data2_pts.append((val2, 100.0*i/(len(data2)-1)))
    data1_pts.append((max_pt, 100))
    data2_pts.append((max_pt, 100))

    # Create the plot
    config = pygal.Config(fill=False, show_legend=True)
    # print("Loading custom CSS file:",CUSTOM_CSS_FILE)
    # config.css.append(CUSTOM_CSS_FILE)
    chart = pygal.XY(config)

    chart.title = "Cumulative Distribution function"  #"CDF: File Size %i -- Record Size %i"%(fsize, rsize)
    chart.x_title = "Throughput"
    chart.y_title = "Percentage of data"
    chart.x_label_rotation = 15
    chart.add("Old 40", data1_pts)
    chart.add("New 40", data2_pts)
    chart.render_to_file(file_name)    

       

# Code for generating the view file for a requested view
if __name__ == "__main__":
    # Get the three primary command line arguments
    group_id = "VirginiaTech"
    project_name = ["Raw_CFQ_NOOP_64_Fread", "F_R_Size_Raw"]
    data_req_base = "[F_size={0}][R_Size={1}][Mode=Fread][Freq=2700000][Throughput=0_&le;_1.0e+20]"
    # All the possible combinations of file and record sizes
    # f_r_sizes = [(64,32),(256,32 ),(1024,32 ),
    #                      (256,128),(1024,128),
    #                                (1024,
    # Ordered in terms of largest to smallest discrepancy variance
    f_r_sizes = [(1024,512), (1024,32), (1024,128),
                 (256,128), (256,32), (64,32)]

    print()
    
    print(HTML_BEGIN)

    plot_pts = []


    # Cycle through the file and record sizes, generating a plot for each
    for f_size, r_size in f_r_sizes:
        # Make the data request string based on the current file and record sizes
        data_req = data_req_base.format(f_size, r_size)

        # Read the data from file to pass to plotting functions,
        # change stdout so that the print statements from
        # "data_for_spec" do not make their way to the terminal
        old_out = sys.stdout
        sys.stdout = open("/dev/null",'w')
        data_old, _,_,_ = data_for_spec(
            group_id, project_name[0], spec_str_to_dict(data_req))
        data_new, _,_,_ = data_for_spec(
            group_id, project_name[1], spec_str_to_dict(data_req))    
        sys.stdout = old_out

        # We only actually want the throughput data from each of these
        study = "Throughput"
        data_old = data_old[study]
        data_new = data_new[study]

        # Process the two sets of 40 numbers to get our nice CDF plot!
        file_name = "Fsize-%i_Rsize-%i.svg"%(f_size, r_size)
        cdf_plot(data_old, data_new, file_name)
        # disparity_plot(data_old, data_new, file_name)

        old_mean, new_mean = data_old.mean(), data_new.mean()
        old_std, new_std = data_old.std(), data_new.std()
        old_var, new_var = data_old.var(), data_new.var()        
        mean_diff = 100.0 * (max(old_mean, new_mean) - min(old_mean, new_mean)) / min(old_mean, new_mean)
        std_diff = 100.0 * (max(old_std, new_std) - min(old_std, new_std)) / min(old_std, new_std)
        var_diff = 100.0 * (max(old_var, new_var) - min(old_var, new_var)) / min(old_var, new_var)

        title = "<strong>File Size %i — Record Size %i</strong>"%(f_size, r_size)
        print(REPORT_ENTRY%(title, mean_diff, std_diff, var_diff, file_name))

        plot_pts.append([f_size, r_size, data_old.var()])
        plot_pts.append([f_size, r_size, data_new.var()])
            

    print()
    print(HTML_MID)
    print(VAR_PLOT % (plot_pts, "File Size", "Record Size", "Variance"))
    print(HTML_END)

# Run command:
# python3 report_code.py > report.html && google-chrome report.html
