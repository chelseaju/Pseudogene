DIR="select_pseudogene_128_v5/LassoValidation_rep54"
SUBDIR=("30X_101L_4A" "30X_101L_6A" "30X_101L_8A"
	"30X_101L_R1A" "30X_101L_R2A" "30X_101L_R3A")


for s in "${SUBDIR[@]}"
do
	X=$s"_reconstructed_x_matrix.txt"
	Y=$s"_validation_y_vector.txt"
	eX=$s"_validation_x_matrix.txt"
	v=$s"_validation"

	R --no-save --slave < lasso_validator_v6.R --args $DIR $X $Y $eX $v

done

