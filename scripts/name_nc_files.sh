
DIR='../data/oceanographic_data/ocean_colour_data/'
for FILE in "$DIR"*.nc
do
	temp_name=$(ncdump $FILE -h | grep :_lastModified | cut -c21-30)
	echo $temp_name$FILE
	without_dir=${FILE/.nc/}
	echo $without_dir$temp_name'.nc'
	mv $FILE $without_dir$temp_name'.nc'
done

