get_er_pt <- function(roc, ...){
  pROC::coords(roc, "best", best.method="closest.topleft")$threshold
}

get_youden_pt <- function(roc, ...){
  pROC::coords(roc, "best", best.method="youden")$threshold
}

get_cz_pt <- function(predicted, actual, pt_seq, confusion_matrices, f_select, ...){
  out <- confusion_matrices[,"Se"]*confusion_matrices[,"Sp"]
  f_select(pt_seq[out==max(out)])
}


get_iu_pt <- function(auc, pt_seq, confusion_matrices, f_select, ...){
  out <- abs(confusion_matrices[, "Se"]-auc) + abs(confusion_matrices[, "Sp"]-auc)
  f_select(pt_seq[out==min(out)])
}


get_nmb_pt <- function(nmb, pt_seq, confusion_matrices, f_select, ...){
  out <- confusion_matrices[,"TN"]*nmb["TN"] + confusion_matrices[,"TP"]*nmb["TP"] + confusion_matrices[,"FN"]*nmb["FN"] + confusion_matrices[,"FP"]*nmb["FP"]
  f_select(pt_seq[out==min(out)])
}
