PRE-REQUISITES:
  python3
  pip3
    scipy
    numpy
    plotly
    colorlover

  In order to meet known working requirements, use:
    $ pip3 install -r requirements.txt
    
USAGE:
  $ python3 compare.py <path to new data CSV file>

  <new_data> must be a path to a comma separated data file 
  with the first row being headers from the following list
  of names:

     Name        Type       Options
     --------    ------     --------
     Machine:    'str'      new
     Store:      'str'      HDD
     Journal:    'str'	    yes
     Hyp:        'str'	    xen
     Hyp_Sched:  'str'	    CFQ, DEAD, NOOP
     VM_Sched:   'str'	    CFQ, DEAD, NOOP
     RCU:        'float'    128
     F_Size:     'float'    64, 256, 1024
     R_Size:     'float'    32, 128, 512
     Threads:    'float'    1, 2, 4, 8, 16, 32, 64, 128, 256
     Mode:       'str'	    Fread, Fwrite, Initialwrite, 
                            Mixedworkload, Pread, Pwrite,
			    Randomread, Randomwrite,
			    Re-read, Read, ReverseRead,
			    Rewrite, Strideread
     Freq:       'float'    1200000, 1400000, 1500000, 1600000, 
                            1800000, 1900000, 2000000, 2100000,
			    2300000, 2400000, 2500000, 2700000,
			    2800000, 2900000, 3000000, 3001000 
     Throughput: 'float'    <float>


  The headers can be in any order, but 'Throughput' *must*
  be included. The values in the columns are assumed to match
  the specified types above. The existing VarSys data will be
  filtered to those rows (data points) which *exactly* match
  the settings for (at least) one of the rows in the provided  
  data file.

  Behavior is not defined for missing values. Excluded columns
  are not used to filter old data (this can be messy, specify
  as many columns as possible for more digestable results).

  Example data provided in 'Resources/sample_data.txt'
  Run example using:
    $ python3 compare.py Resources/sample_data.txt

