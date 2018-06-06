mkdir /ahg/regevdata/projects/MerckIT/GitHub/
cd /ahg/regevdata/projects/MerckIT/GitHub/

mkdir Code
mkdir Data
mkdir Output
mkdir Results

cd /ahg/regevdata/projects/MerckIT/GitHub/Results
mkdir CellTypes
mkdir Predictors
mkdir Resistance
mkdir RevisionResults
mkdir Signatures
mkdir Spatial
mkdir Resistance/Exclusion
mkdir Resistance/Functional
mkdir Resistance/PostTreatment

cp /ahg/regevdata/projects/MerckIT/GitHubPrivate/Results/Resistance/PostTreatment/post.treatment.subsample.de.RData Resistance/PostTreatment/

cd /ahg/regevdata/projects/MerckIT/GitHub/Data/
mkdir scData
mkdir PublicData
mkdir ValidationCohorts

mkdir /ahg/regevdata/projects/MerckIT/GitHub/Output/Tables
mkdir /ahg/regevdata/projects/MerckIT/GitHub/Output/Figures

F=(Mel.all.data.QC.rds Mel.all.data.QC_genes.rds Mel.malignant.rds Mel.malignant_genes.rds Mel.T.CD4.QC.rds Mel.T.CD8.QC.rds)

for f in ${F[*]}
do
echo "copying $f"
cp /ahg/regevdata/projects/MerckIT/GitHubPrivate/Data/scData/$f /ahg/regevdata/projects/MerckIT/GitHub/Data/scData/
done

cp /ahg/regevdata/projects/MerckIT/GitHubPrivate/Data/PublicData/T.cell.state.markers.rds /ahg/regevdata/projects/MerckIT/GitHub/Data/PublicData/
cp /ahg/regevdata/projects/MerckIT/GitHubPrivate/Data/PublicData/mouse_human_mapping.RData /ahg/regevdata/projects/MerckIT/GitHub/Data/PublicData/
cp /ahg/regevdata/projects/MerckIT/GitHubPrivate/Data/PublicData/Sanger.Garnett.data.rds /ahg/regevdata/projects/MerckIT/GitHub/Data/PublicData/