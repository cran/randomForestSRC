\name{breast}
\docType{data}
\alias{breast}
\title{Wisconsin Prognostic Breast Cancer Data}
\description{
  Recurrence of breast cancer from 198 breast cancer patients, all of
  which exhibited no evidence of distant metastases at the time of
  diagnosis.  The first 30 features of the data describe characteristics
  of the cell nuclei present in the digitized image of a fine needle
  aspirate (FNA) of the breast mass.
}
\format{A data frame containing:
  \tabular{ll}{

    status \tab factor with levels N (nonrecurrent) and R (recurrent)
    indicating the patients outcome\cr  
    mean_radius \tab radius (mean of distances from center to points on the perimeter) (mean)\cr
    mean_texture \tab texture (standard deviation of gray-scale values) (mean)\cr
    mean_perimeter \tab perimeter (mean)\cr
    mean_area \tab area (mean)\cr
    mean_smoothness \tab smoothness (local variation in radius lengths) (mean)\cr
    mean_compactness \tab compactness (mean)\cr
    mean_concavity \tab concavity (severity of concave portions of the contour) (mean)\cr
    mean_concavepoints \tab concave points (number of concave portions of the contour) (mean)\cr
    mean_symmetry \tab symmetry (mean)\cr
    mean_fractaldim \tab fractal dimension (mean)\cr
    SE_radius \tab radius (mean of distances from center to points on the perimeter) (SE)\cr
    SE_texture \tab texture (standard deviation of gray-scale values) (SE)\cr
    SE_perimeter \tab perimeter (SE)\cr
    SE_area \tab area (SE)\cr
    SE_smoothness \tab smoothness (local variation in radius lengths) (SE)\cr
    SE_compactness \tab compactness (SE)\cr
    SE_concavity \tab concavity (severity of concave portions of the contour) (SE)\cr
    SE_concavepoints \tab concave points (number of concave portions of the contour) (SE)\cr
    SE_symmetry \tab symmetry (SE)\cr
    SE_fractaldim \tab fractal dimension (SE)\cr
    worst_radius \tab radius (mean of distances from center to points on the perimeter) (worst)\cr
    worst_texture \tab texture (standard deviation of gray-scale values) (worst)\cr
    worst_perimeter \tab perimeter (worst)\cr
    worst_area \tab area (worst)\cr
    worst_smoothness \tab smoothness (local variation in radius lengths) (worst)\cr
    worst_compactness \tab compactness  (worst)\cr
    worst_concavity \tab concavity (severity of concave portions of the contour) (worst)\cr
    worst_concavepoints \tab concave points (number of concave portions of the contour) (worst)\cr
    worst_symmetry \tab symmetry (worst)\cr
    worst_fractaldim \tab fractal dimension (worst)\cr
    tsize \tab diameter of the excised tumor in centimeters\cr
    pnodes \tab number of positive axillary lymph nodes observed at time of surgery\cr
  }
}
\source{
  The data were obtained from the UCI machine learning repository, see
  \url{http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Prognostic)}.
}
\examples{
\dontrun{
data(breast, package = "randomForestSRC")
breast.obj <- rfsrc(status ~ ., data = breast, nsplit = 10)
print(breast.obj)
plot(breast.obj)
}}
\keyword{datasets}
