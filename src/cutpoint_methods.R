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

fx_total_nmb <- function(tn, tp, fn, fp, utility_tp, utility_tn, cost_fp, cost_fn, ...){
  total_nmb <- utility_tp * tp + utility_tn * tn + cost_fp * fp + cost_fn * fn
  total_nmb <- matrix(total_nmb, ncol = 1)
  colnames(total_nmb) <- "total_nmb"
  total_nmb
}

fx_prod_sn_sp <- function(tp, fp, tn, fn, ...) {
  sesp <- cutpointr:::sens_spec(tp, fp, tn, fn)
  prod_sn_sp <- sesp[,1]*sesp[,2]
  prod_sn_sp <- matrix(prod_sn_sp, ncol=1)
  colnames(prod_sn_sp) <- "prod_sn_sp"
  prod_sn_sp
}
