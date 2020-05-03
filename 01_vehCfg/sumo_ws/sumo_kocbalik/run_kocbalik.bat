sumo -c osm.sumocfg --device.fcd.explicit veh8,veh12,veh13
move ".\kocbalik.xml" ".\xml_convert"
cd xml_convert
python xml2csv.py -s "," kocbalik.xml
move ".\kocbalik.csv" "..\data\kocbalik.csv"
del ".\kocbalik.xml"
pause