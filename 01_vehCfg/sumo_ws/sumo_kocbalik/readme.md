The "run_kocbalik.bat" file was derived from the default "run.bat" file, which originally just instantiated the GUI, for this purpose. "run_kocbalik.bat" instantiates the command line version of the SUMO package and on top of running the SUMO simulation, it uses Windows commands to copy the generated [FCD](https://sumo.dlr.de/docs/Simulation/Output/FCDOutput.html) XML file to an xml_convert directory and calls the xml2csv.py Python tool to convert it to CSV. 

The simulation step time (i.e., vehicle simulation sampling period in our books), and the xml file output name ("kocbalik.xml" in our case), are important parameters and are edited manually inside the osm.sumocfg file, which is actually in XML format (so just open it in your text editor of choice)

veh8, veh12 and veh13 were manually selected after observing a GUI run of this SUMO scenario for our simulation purposes, so there are only FCD CSV files for those three. Two of those vehicles were later selected as ego and target vehicles for a simulation provided in the article under review, as shown below

<img src="../../../99_doc/sumo_kocbalik.png" alt="Drawing"/>