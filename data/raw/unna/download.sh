for year in {1960..2010} 
do
	curl "http://data.un.org/Handlers/DownloadHandler.ashx?DataFilter=group_code:201;fiscal_year:$year&DataMartId=SNA&Format=csv" > temp.zip
    unzip -p temp.zip > $year.csv
done
rm temp.zip
