#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;


IntegerVector pure_erase(IntegerVector x, int ind) {
  IntegerVector z=clone(x);
  z.erase(ind);
  return(z);
}

NumericVector pure_erase(NumericVector x, int ind) {
  NumericVector z=clone(x);
  z.erase(ind);
  return(z);
}

NumericMatrix mmultcpp( NumericMatrix m1, NumericMatrix m2){
  if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
  NumericMatrix out(m1.nrow(),m2.ncol());
  NumericVector rm1, cm2;
  for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
  return out;
}

NumericVector mmultcpp(NumericMatrix mm,NumericVector vv){
  NumericVector out(mm.nrow());
  NumericVector rm1;
  for (int i=0; i<mm.nrow();i++) {
    rm1=mm(i,_);
    out[i]=std::inner_product(rm1.begin(),rm1.end(),vv.begin(),0.0);
  }
  return out;
}


double sumveccpp(NumericVector v1) {
  int n=v1.size();
  double res=0;
  for(int i=0;i<n;i++) {
    res += v1[i];
  }
  return res;
}

NumericVector combine(const List& list)
{
  std::size_t n = list.size();
  
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  
  // Allocate the vector
  NumericVector output = no_init(total_length);
  
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    
    // Update the index
    index += el.size();
  }
  
  return output;
  
}


// NumericVector beta_mucpp(NumericVector t_s) {
//   /* double a_mu=3;
//    double b_mu=5;
//    double c_mu=1;
//    NumericVector rt=(a_mu*(1-b_mu/(5+(c_mu*t_s))));*/
//   NumericVector rt=20.0-(4.0/(1.0+(0.2*t_s)));
//   return rt;
// }

NumericMatrix beta_piecewise(NumericVector t_s, NumericVector knots, double t_end, NumericVector beta) {
  int K=beta.size()-1;
  int L=t_s.size();
  double time_int;
  double t_temp;
  int k_start;
  double sum_intercept;
  double res_ti;
  double slope;
  NumericVector res_vec(L);
  NumericVector intercept_vec(L);
  NumericVector slope_vec(L);
  NumericVector knots_start = knots[seq(0,K-1)];
  NumericVector knots_end = knots[seq(1,K)];
  NumericVector knots_start_vec(L);
  // IntegerVector k_index=seq(0,K-1);
  // NumericVector k_numeric=as<NumericVector>(k_index);
  // NumericVector knot_start=k_numeric*time_int;
  // NumericVector knot_end=time_int * (k_numeric+1);
  IntegerVector k_start_v;
  IntegerVector index=seq(0,K-1);
  // Rcout << " t_s " << t_s << "\n";
  for (int i=0;i<L;i++) {
    t_temp=t_s[i];
    // Rcout << " t_temp " << t_temp << "\n";
    if (t_temp>=t_end) {
      t_temp=t_end;
      k_start=K-1;
    } else if (t_temp<=knots[0]) {
      t_temp = knots[0];
      k_start=0;
    } else {
      k_start_v=index[t_temp>=knots_start & t_temp<=knots_end];
      // Rcout << " k_start_v " << k_start_v << "\n";
      k_start=k_start_v[0];
    }
    // Rcout << " k_start " << k_start << " i " << i << " L " << L << "\n";
    
    
    // double t_start=t_start_v[0];
    sum_intercept=beta[0];
    for (int j=0;j<k_start;j++) {
      time_int=knots[j+1]-knots[j];
      sum_intercept=sum_intercept+beta[j+1]*time_int;
    }
    // Rcout << " sum_intercpt " << sum_intercept << " knots[k_start] " << knots[k_start]  << " beta[k_start+1] " << beta[k_start+1] << "\n";
    res_ti=sum_intercept + beta[k_start+1] * (t_temp-knots[k_start]);
    res_vec[i]=res_ti;
    slope=beta[k_start+1];
    if (t_temp>=t_end) {
      sum_intercept=res_ti;
      slope=0;
    }
    intercept_vec[i]=sum_intercept;
    slope_vec[i]=slope;
    knots_start_vec[i]=knots[k_start];
    // Rcout << " res_ti " << res_ti << "\n";
    // Rcout << " time_int " << time_int  << " sum_intercept " << sum_intercept << " k_start " << k_start << " beta_used " << beta[k_start+1] << "\n"; 
  }
  NumericMatrix res_mat(L,4);
  res_mat(_,0)=res_vec;
  res_mat(_,1)=intercept_vec;
  res_mat(_,2)=slope_vec;
  res_mat(_,3)=knots_start_vec;
  return res_mat;
}


// NumericVector beta_Acpp(NumericVector t_s) {
//   /* double a_mu=3;
//    double b_mu=5;
//    double c_mu=1;
//    NumericVector rt=(a_mu*(1-b_mu/(5+(c_mu*t_s))));*/
//   NumericVector rt=4.0*exp(-0.23*t_s-0.92);
//   return rt;
// }

// double alphaAcpp (double x) {
//   return ( -x/12 + 2.5);
// }



NumericVector muicpp(NumericVector tstari, double tau, double betamu, double betaA, NumericVector PsiX, NumericVector Xi, double Ai,
                     double Ui, NumericVector beta_mu_vec, NumericVector beta_A_vec, 
                     NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, NumericVector knots, double t_end) {
  
  double t_end_1=max(knots);
  NumericVector Xcol1=beta_piecewise(tstari,knots,t_end_1,beta_mu_vec)(_,0);
  NumericVector Xcol2=Ai*beta_piecewise(tstari,knots,t_end_1,beta_A_vec)(_,0);
  NumericVector Xcol3=beta_piecewise(tstari,knots,t_end_1,beta_mu_Ui_vec)(_,0);
  NumericVector Xcol4=Ai*beta_piecewise(tstari,knots,t_end_1,beta_A_Ui_vec)(_,0);
  
  NumericVector m = betamu*Xcol1 + betaA*Xcol2+ (Xcol3+Xcol4)*Ui + sum(PsiX*Xi);
  // Rcout << " beta* " << mmultcpp(XitstarT,beta) << " betamu_int " <<beta_mu_int <<  " betamu_slope " << beta_mu_slope <<" betaA_int " <<beta_A_int <<  " betaA_slope " << beta_A_slope << "\n";
  return m;
}




NumericVector int1cpp (double ti,double gamma0,double gamma1,NumericVector gammaX,double gammaU,NumericVector Xi,double Ai,double Ui,double sty_end,
                       NumericVector drop_knots, NumericVector drop_par_0, NumericVector drop_par_1){
  NumericVector drop_knots_start=clone(drop_knots);
  drop_knots_start.push_front(0);
  NumericVector drop_knots_end=clone(drop_knots);
  drop_knots_end.push_back(sty_end);
  
  IntegerVector find_surv_loc = ifelse((ti > drop_knots_start) & (ti <= drop_knots_end),1,0);
  int surv_loc=which_max(find_surv_loc);
  if (sum(find_surv_loc)==0) {
    surv_loc=drop_knots.size();
  }
  
  drop_par_0.push_front(gamma0);
  drop_par_1.push_front(gamma1);
  double int1_base=0;
  for (int i=0;i<surv_loc;i++) {
    int1_base=int1_base + (drop_knots_start[i+1] - drop_knots_start[i])*exp((gammaU*Ui)+drop_par_0[i]+(Ai*drop_par_1[i])+sum(Xi*gammaX));
  }
  double rtn = int1_base + exp((gammaU*Ui)+drop_par_0[surv_loc]+(Ai*drop_par_1[surv_loc])+sum(Xi*gammaX))*(ti-drop_knots_start[surv_loc]);
  double base = int1_base - exp((gammaU*Ui)+drop_par_0[surv_loc]+(Ai*drop_par_1[surv_loc])+sum(Xi*gammaX)) * drop_knots_start[surv_loc];
  NumericVector res={rtn,drop_par_0[surv_loc],drop_par_1[surv_loc],base};
  return(res);
}

NumericVector int2cpp (double ti,double alpha0,double alpha1,NumericVector alphaX,double alphaU,NumericVector Xi,double Ai,double Ui,double sty_end,
                       NumericVector surv_knots, NumericVector surv_par_0, NumericVector surv_par_1){
  
  NumericVector surv_knots_start=clone(surv_knots);
  surv_knots_start.push_front(0);
  NumericVector surv_knots_end=clone(surv_knots);
  surv_knots_end.push_back(sty_end);
  
  IntegerVector find_surv_loc = ifelse((ti > surv_knots_start) & (ti <= surv_knots_end),1,0);
  int surv_loc=which_max(find_surv_loc);
  if (sum(find_surv_loc)==0) {
    surv_loc=surv_knots.size();
  }
  surv_par_0.push_front(alpha0);
  surv_par_1.push_front(alpha1);
  double int1_base=0;
  for (int i=0;i<surv_loc;i++) {
    int1_base=int1_base + (surv_knots_start[i+1] - surv_knots_start[i])*exp((alphaU*Ui)+surv_par_0[i]+(Ai*surv_par_1[i])+sum(Xi*alphaX));
  }
  double rtn = int1_base + exp((alphaU*Ui)+surv_par_0[surv_loc]+(Ai*surv_par_1[surv_loc])+sum(Xi*alphaX))*(ti-surv_knots_start[surv_loc]);
  double base = int1_base - exp((alphaU*Ui)+surv_par_0[surv_loc]+(Ai*surv_par_1[surv_loc])+sum(Xi*alphaX)) * surv_knots_start[surv_loc];
  
  NumericVector res={rtn,surv_par_0[surv_loc],surv_par_1[surv_loc],base};
  return(res);
}

// 
// NumericMatrix find_max(double (*f)(double ,List ), List par_list, double left_start, double right_start) {
//   int num_interv=100;
//   double interv=20.0/num_interv;
//   double left_inp= left_start;
//   double left_point_start = left_start;
//   double left_value;
//   double new_n;
//   // double new_n1;
//   double left;
//   double right;
//   double frst;
//   double M_temp;
//   double M;
//   double xa;
//   double xb;
//   // double M1;
//   NumericVector temp_save(num_interv);
//   for (int p=0;p<num_interv;p++) {
//     left_value= (*f)(left_inp, par_list);
//     left_inp=left_inp+interv;
//     temp_save[p]=left_value;
//   }
//   // int max_ind=which_max(temp_save);
//   double max_temp=max(temp_save);
//   //  left= max_ind*interv;
//   //  right = (max_ind+1)*interv;
//   IntegerVector index_seq=seq(0,99);
//   
//   NumericVector pot_cluster_start;
//   NumericVector pot_cluster_end;
//   //temp_save len 100
//   index_seq=index_seq[temp_save>max_temp-20];
//   
//   temp_save=temp_save[temp_save>max_temp-20];
//   // NumericVector temp_save_ori=clone(temp_save);
//   // temp_save_ori.erase(temp_save.size()-1);
//   NumericVector temp_save_start=pure_erase(temp_save,temp_save.size()-1);
//   NumericVector temp_save_end=pure_erase(temp_save,0);
//   
//   IntegerVector index_start=pure_erase(index_seq,index_seq.size()-1);
//   IntegerVector index_end=pure_erase(index_seq,0);
//   LogicalVector pot_cluster_indi=(index_end-index_start)>1;
//   
//   IntegerVector index_start_with=index_end[pot_cluster_indi];
//   IntegerVector index_end_with=index_start[pot_cluster_indi];
//   
//   index_start_with.push_front(index_seq[0]);
//   index_end_with.push_back(index_seq[index_seq.size()-1]);
//   
//   NumericVector M_vec(index_start_with.size());
//   NumericVector xa_vec(index_start_with.size());
//   NumericVector xb_vec(index_start_with.size());
//   
//   // Rcout << " temp_save " << temp_save << "\n"; 
//   // Rcout << " index_seq " << index_seq << "\n";
//   // Rcout << " index_start_with " << index_start_with << "\n";
//   // Rcout << " index_end_with " << index_end_with << "\n";
//   double left_cluster;
//   double right_cluster;
//   for (int i=0;i<index_start_with.size();i++) {
//     left_cluster=left_point_start + index_start_with[i]*interv;
//     right_cluster=left_point_start + index_end_with[i]*interv;
//     left_cluster = left_cluster-interv;
//     // if (left_cluster<0) {
//     //   left_cluster=0;
//     // }
//     right_cluster = right_cluster+interv;
//     
//     left=left_cluster;
//     right=right_cluster;
//     // new_n=left;
//     frst=1;
//     // Rcout << " left_inp " << left_inp << " max_ind " << max_ind << " max_temp " << max(temp_save) << " temp_save " << temp_save << "\n";
//     while (std::abs(frst)>0.001 && std::abs(left-right)>1e-6) {
//       new_n=(left+right)/2;
//       frst=((*f)(new_n+1e-3, par_list) - (*f)(new_n, par_list))/1e-3;
//       // Rcout << " new_n_loop " << new_n << "\n";
//       if (frst>0) {
//         left=new_n;
//       } else {
//         right=new_n;
//       }
//     }
//     
//     M = (*f)(new_n, par_list);
//     // M1 = li1cpp(ti,yi,tau,sigma,PsiX,gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, betamu,
//     //             betaA, Xi, Ai, tstari,beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , 0);
//     // if (M1>M) {
//     //   M=M1;
//     //   new_n=0;
//     // }
//     // Rcout << " M " << M << " new_n " << new_n << "\n";
//     left=new_n;
//     right=right_cluster;
//     
//     M_temp=(*f)(right, par_list);
//     
//     if (M_temp>M-20) {
//       xb=right;
//     } else {
//       M_temp=M+1;
//       while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>1e-6) {
//         xb=(left+right)/2;
//         M_temp=(*f)(xb, par_list);
//         
//         if (M_temp>M-20) {
//           left=xb;
//         } else {
//           right=xb;
//         }
//       }
//     }
//     left=left_cluster;
//     right=new_n;
//     M_temp=(*f)(left, par_list);
//     
//     if (M_temp>M-20) {
//       xa=left;
//     } else {
//       M_temp=M+1;
//       while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>1e-6) {
//         xa=(left+right)/2;
//         M_temp=(*f)(xa, par_list);
//         if (M_temp>M-20) {
//           right=xa;
//         } else {
//           left=xa;
//         }
//       }
//     }
//     M_vec[i]=M;
//     xa_vec[i]=xa;
//     xb_vec[i]=xb;
//   }
//   NumericMatrix result(index_start_with.size(),3);
//   result (_,0) = M_vec;
//   result (_,1) = xa_vec;
//   result (_,2) = xb_vec;
//   // Rcout << " result " << result << "\n";
//   return result;
// }


NumericMatrix find_max(double (*f)(double ,List ), List par_list, double left_start, double right_start) {
  
  int num_interv=100;
  double right_limit=right_start;
  double left_inp=left_start;
  double interv=(right_limit - left_inp)/num_interv;
  double left_start_point = left_start;
  
  double left_value;
  double new_n;
  double new_n1=0;
  double left;
  double right;
  double frst;
  double M_temp;
  double M;
  double M1;
  double xa;
  double xb;
  NumericVector temp_save_temp(num_interv);
  // NumericVector temp_left_temp(num_interv);
  List temp_save_list;
  // List temp_left_list;
  // Rcout << " find2_start1 " << "\n";
  
  int max_ind=num_interv-1;
  while (max_ind==num_interv-1) {
    for (int p=0;p<num_interv;p++) {
      left_value=(*f)(left_inp, par_list);
      // Rcout << " left_inp " << left_inp << " left_value " << left_value << "\n";
      temp_save_temp[p]=left_value;
      // temp_left_temp[p]=left_inp;
      left_inp=left_inp+interv;
    }
    // right_limit=right_limit*2;
    // interv=right_limit/num_interv;
    max_ind=which_max(temp_save_temp);
    temp_save_list.push_back(clone(temp_save_temp));
    // temp_left_list.push_back(temp_left_temp);
  }
  
  NumericVector temp_save1=combine(temp_save_list);
  max_ind = which_max(temp_save1);
  
  while (max_ind==0) {
    left_start_point = left_start_point - (right_start - left_start);
    left_inp = left_start_point;
    for (int p=0;p<num_interv;p++) {
      left_value=(*f)(left_inp, par_list);
      // Rcout << " left_inp " << left_inp << " left_value " << left_value << "\n";
      temp_save_temp[p]=left_value;
      // temp_left_temp[p]=left_inp;
      left_inp=left_inp+interv;
    }
    // right_limit=right_limit*2;
    // interv=right_limit/num_interv;
    max_ind=which_max(temp_save_temp);
    temp_save_list.push_front(clone(temp_save_temp));
    // temp_left_list.push_back(temp_left_temp);
  }
  
  
  NumericVector temp_save=combine(temp_save_list);
  // NumericVector temp_left=combine(temp_left_list);
  
  
  // Rcout << " find2_start2 " << "\n";
  
  double max_temp=max(temp_save);
  max_ind = which_max(temp_save);
  //  left= max_ind*interv;
  //  right = (max_ind+1)*interv;
  IntegerVector index_seq=seq(0,temp_save.size()-1);
  IntegerVector index_seq2=seq(0,temp_save.size()-1);
  
  // Rcout << " temp_save " << temp_save << " max_temp " << max_temp << "\n";
  
  NumericVector pot_cluster_start;
  NumericVector pot_cluster_end;
  //temp_save len 100
  // Rcout << " index_seq " << index_seq << "\n";
  
  index_seq=index_seq[temp_save>max_temp-20];
  
  // Rcout << " index_seq " << index_seq << "\n";
  IntegerVector index_start=pure_erase(index_seq,index_seq.size()-1);
  IntegerVector index_end=pure_erase(index_seq,0);
  LogicalVector pot_cluster_indi=(index_end-index_start)>1;
  
  // Rcout << " index_start " << index_start  << " index_end " << index_end << " pot_cluster_indi " << pot_cluster_indi << "\n";
  
  IntegerVector index_start_with=index_end[pot_cluster_indi];
  IntegerVector index_end_with=index_start[pot_cluster_indi];
  // Rcout << " index_start_with " << index_start_with << " index_end_with " << index_end_with << "\n";
  
  
  index_start_with.push_front(index_seq[0]);
  index_end_with.push_back(index_seq[index_seq.size()-1]);
  // Rcout << " index_start_with " << index_start_with << " index_end_with " << index_end_with << "\n";
  
  NumericVector M_vec(index_start_with.size());
  NumericVector xa_vec(index_start_with.size());
  NumericVector xb_vec(index_start_with.size());
  // Rcout << " find2_start3 " << "\n";
  
  // Rcout << " temp_save " << temp_save << "\n";
  // Rcout << " index_seq " << index_seq << "\n";
  // Rcout << " index_start_with " << index_start_with << "\n";
  // Rcout << " index_end_with " << index_end_with << "\n";
  double left_cluster;
  double right_cluster;
  NumericVector temp_save2;
  IntegerVector temp_index2;
  int sub_max_indi;
  int max_indi_temp;
  for (int i=0;i<index_start_with.size();i++) {
    
    temp_save2=temp_save[Rcpp::Range(index_start_with[i],index_end_with[i])];
    temp_index2=index_seq2[Rcpp::Range(index_start_with[i],index_end_with[i])];
    sub_max_indi=which_max(temp_save2);
    max_indi_temp=temp_index2[sub_max_indi];
    
    
    left_cluster= left_start_point + max_indi_temp*interv;
    right_cluster=left_cluster+interv;
    // left_cluster = left_cluster-interv;
    // if (left_cluster<0) {
    //   left_cluster=0;
    // }
    // right_cluster = right_cluster+interv;
    
    left=left_cluster;
    right=right_cluster;
    // new_n=left;
    frst=1;
    // Rcout << " left_inp " << left_inp << " max_ind " << max_ind << " max_temp " << max(temp_save) << " temp_save " << temp_save << "\n";
    // Rcout << " left_cluster " << left_cluster << " right_cluster " << right_cluster << "\n";
    while (std::abs(frst)>0.001 && std::abs(left-right)>1e-6) {
      new_n=(left+right)/2;
      frst=((*f)(new_n+1e-6, par_list) - (*f)(new_n, par_list))/1e-6;
      // Rcout << " new_n_loop " << new_n << " left_loop " << left << " right_loop " << right << " frst " << frst  << "\n";
      if (frst>0) {
        left=new_n;
      } else {
        right=new_n;
      }
    }
    
    
    M = (*f)(new_n, par_list);
    // Rcout << " find2_start4 " << "\n";
    
    // Rcout << " new_n " << new_n << " M " << M << "\n";
    right_cluster=left_start_point + max_indi_temp*interv;
    left_cluster=right_cluster-interv;
    // left_cluster = left_cluster-interv;
    // if (left_cluster<0) {
    //   left_cluster=0;
    // }
    // right_cluster = right_cluster+interv;
    
    left=left_cluster;
    right=right_cluster;
    // new_n=left;
    frst=1;
    new_n1=new_n;
    // Rcout << " left_inp " << left_inp << " max_ind " << max_ind << " max_temp " << max(temp_save) << " temp_save " << temp_save << "\n";
    while (std::abs(frst)>0.001 && std::abs(left-right)>1e-6) {
      new_n1=(left+right)/2;
      frst=((*f)(new_n1+1e-6, par_list) - (*f)(new_n1, par_list))/1e-6;
      // Rcout << " new_n_loop " << new_n1 << " left_loop " << left << " right_loop " << right << " frst " << frst  << "\n";
      if (frst>0) {
        left=new_n1;
      } else {
        right=new_n1;
      }
    }
    
    
    M1 = (*f)(new_n1, par_list);
    // Rcout << " find2_start5 " << "\n";
    
    if (M1>M) {
      M=M1;
      new_n=new_n1;
    }
    
    // M1 = li2cpp(ti,yi,tau,sigma,PsiX,gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, betamu,
    //             betaA, Xi, Ai, tstari, 0,longi_t,beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ,tol,nodes,weights);
    // 
    // if (M1>M) {
    //   M=M1;
    //   new_n=0;
    // }
    
    // Rcout << " M " << M << " new_n " << new_n << "\n";
    left_cluster=left_start_point + index_start_with[i]*interv - interv;
    right_cluster=left_start_point + index_end_with[i]*interv + interv;
    // if (left_cluster<0) {
    //   left_cluster=0;
    // }
    left=new_n;
    right=right_cluster;
    
    M_temp=(*f)(right, par_list);
    
    if (M_temp>M-20) {
      xb=right;
    } else {
      M_temp=M+1;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>1e-6) {
        xb=(left+right)/2;
        M_temp=(*f)(xb, par_list);
        if (M_temp>M-20) {
          left=xb;
        } else {
          right=xb;
        }
      }
    }
    // Rcout << " find2_start6 " << "\n";
    
    left=left_cluster;
    right=new_n;
    M_temp=(*f)(left, par_list);
    
    if(M_temp>M-20) {
      xa=left;
    } else {
      M_temp=M+1;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>1e-6) {
        xa=(left+right)/2;
        M_temp=(*f)(xa, par_list);
        
        if (M_temp>M-20) {
          right=xa;
        } else {
          left=xa;
        }
      }
    }
    M_vec[i]=M;
    xa_vec[i]=xa;
    xb_vec[i]=xb;
  }
  // Rcout << " find2_end " << "\n";
  NumericMatrix result(index_start_with.size(),3);
  result (_,0) = M_vec;
  result (_,1) = xa_vec;
  result (_,2) = xb_vec;
  // Rcout << " result " << result << "\n";
  
  
  
  return result;
}


double ffdenfunc(double y, double (*f)(double, List), List par_list, double M, double xa, double xb) {
  double u = -1; 
  double v = 1 ; 
  double duv = (xb-xa)/(v-u);
  return duv * exp(-M + (*f)(xa + duv*(y-u), par_list));
}


double quadinfIntcppdenfunc(double tol, List nodes, List weights, double (*f)(double, List), List par_list, double M,double xa, double xb) {
  double Q;
  double h=0.5;
  NumericVector x=nodes[0];
  NumericVector w=weights[0];
  double s= w[6] * ffdenfunc(x[6], (*f) ,par_list, M,xa,xb);
  for (int j=0;j<6;j++) {
    s=s + w[j] * (ffdenfunc(x[j], (*f) ,par_list, M,xa,xb) + ffdenfunc(-x[j], (*f) ,par_list ,M,xa,xb));
  }
  Q=s*h;
  //int k1;
  double delta1;
  double newQ;
  for (int k=1;k<7;k++) {
    x=nodes[k];
    w=weights[k];
    s=0;
    for (int j=0;j<w.size();j++) {
      s=s+w[j]*(ffdenfunc(x[j], (*f) ,par_list, M,xa,xb) + ffdenfunc(-x[j], (*f) ,par_list, M,xa,xb));
    }
    h = h/2;
    newQ=s*h + Q/2;
    delta1 = std::abs(newQ-Q);
    Q=newQ;
    if (delta1<tol) {
      break;
    }
  }
  return Q;
}



double li1cpp(double ti,NumericVector yi,double tau, NumericVector PsiX,double gamma0,double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
              double betaA,NumericVector Xi,double Ai,NumericVector tstari,
              NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
              NumericVector knots, double t_end, double sty_end, 
              NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double Ui) {
  int ni=tstari.size();
  double sigma=1;
  double trm1 = -ni*log(tau);
  
  double Ui_a=Ui;
  NumericVector mu_i=muicpp(tstari,tau,betamu,betaA,PsiX,Xi,Ai,Ui_a,beta_mu_vec,beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, knots, t_end);
  double trm2= -sumveccpp((yi-mu_i)*(yi-mu_i)) / (2*tau*tau);
  NumericVector alpha_vec=int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1);
  double trm3= - alpha_vec[0];
  double trm4= alpha_vec[1];
  double trm5= alpha_vec[2]*Ai;
  double trm6=Ui_a*(alphaU-Ui_a/(2*sigma*sigma));
  //  double trm61= -Ui_a*Ui_a/(2*sigma*sigma);
  double trm7=sum(Xi*alphaX);
  double trm8= - int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1)[0];
  double trm9= -log(sigma);
  double terms=trm1+trm2+trm3+trm4+trm5+trm6+trm7+trm8+trm9;
  
  //  double terms=trm1+trm2+trm61+trm9;
  return terms;
}


double li1cpp_shell(double Ui, List par_list) {
  double ti=par_list[0];
  NumericVector yi=par_list[1];
  NumericVector thetam=par_list[2];
  NumericVector PsiX=par_list[3];
  NumericVector Xi=par_list[4];
  double Ai=par_list[5];
  NumericVector tstari=par_list[6];
  NumericVector beta_mu_vec=par_list[7];
  NumericVector beta_A_vec=par_list[8];
  NumericVector knots=par_list[9];
  double t_end=par_list[10];
  double sty_end=par_list[11];
  NumericVector surv_knots=par_list[12];
  NumericVector drop_knots=par_list[13];
  NumericVector surv_par_0=par_list[14];
  NumericVector surv_par_1=par_list[15];
  NumericVector drop_par_0=par_list[16];
  NumericVector drop_par_1=par_list[17];
  
  List add_elements=par_list[18];
  
  NumericVector alphaX=add_elements[0];
  NumericVector gammaX =add_elements[1];
  NumericVector beta_mu_Ui_vec=add_elements[2];
  NumericVector beta_A_Ui_vec=add_elements[3];
  
  double betamu=thetam[0];
  double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];

  double alphaU=thetam[4];
  double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];

  double gammaU=thetam[8];

  
  double res=li1cpp(ti,yi,tau,PsiX,gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, betamu,
                    betaA, Xi, Ai, tstari,beta_mu_vec, beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec,
                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , Ui);
  return res;
  
}

double trm3outcpp(double s, double ti,NumericVector yi,double tau,NumericVector PsiX,double gamma0,double gamma1,
                  NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
                  double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
                  NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                  NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1) {
  NumericVector tstar_s = s-longi_t;
  // Rcout << " s " << s << " Ui " << Ui << "\n";
  double Ui_a = Ui;
  NumericVector mu_i=muicpp(tstar_s,tau,betamu,betaA,PsiX,Xi,Ai,Ui_a,beta_mu_vec,beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end);
  double trm_inner31= -sumveccpp((yi-mu_i)*(yi-mu_i))/(2*tau*tau);
  // double trm_inner31= -((yi-mu_i)*(yi-mu_i))[12]/(2*tau*tau);
  NumericVector trm32_3_vec=int2cpp(s,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1);
  double trm_inner32= trm32_3_vec[1]+(trm32_3_vec[2]*Ai);
  double trm_inner33= -(trm32_3_vec[0] - int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1)[0]);
  double rtn=trm_inner31+trm_inner32+trm_inner33;
  // Rcout << " inneroutp " << rtn << " ui " << Ui << "\n";
  // Rcout << " mui " << mu_i << " trm_inner31 " << trm_inner31 << " trm_inner33 " << trm_inner33 << "\n";
  return rtn;
  
}

NumericVector intergrate_trm3outcpp(double xa, double xb, double ti,NumericVector yi,double tau,double sigma,NumericVector PsiX,double gamma0,double gamma1,
                                    NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
                                    double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
                                    NumericVector beta_mu_vec, NumericVector beta_A_vec,NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                                    NumericVector knots, double t_end, double sty_end, 
                                    NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1) {
  
  NumericVector tstar_s=(xb+xa)/2 -longi_t;
  double t_end_1 = max(knots);
  NumericMatrix inter_slope_mu=beta_piecewise(tstar_s,knots,t_end_1,beta_mu_vec);
  NumericMatrix inter_slope_A=beta_piecewise(tstar_s,knots,t_end_1,beta_A_vec);
  NumericMatrix inter_slope_mu_Ui=beta_piecewise(tstar_s,knots,t_end_1,beta_mu_Ui_vec);
  NumericMatrix inter_slope_A_Ui=beta_piecewise(tstar_s,knots,t_end_1,beta_A_Ui_vec);
  
  NumericVector inter_mu=betamu*inter_slope_mu(_,1);
  NumericVector slope_mu=betamu*inter_slope_mu(_,2);
  NumericVector knots_start=inter_slope_mu(_,3);
  NumericVector inter_A=betaA*inter_slope_A(_,1);
  NumericVector slope_A=betaA*inter_slope_A(_,2);
  
  NumericVector inter_mu_Ui=inter_slope_mu_Ui(_,1);
  NumericVector slope_mu_Ui=inter_slope_mu_Ui(_,2);
  
  NumericVector inter_A_Ui=inter_slope_A_Ui(_,1);
  NumericVector slope_A_Ui=inter_slope_A_Ui(_,2);
  
  // NumericVector tstar_tend_comp = ifelse( tstar_s>=t_end, 1.0, 0.0);
  NumericVector beta_0=((inter_mu+inter_mu_Ui*Ui) - (slope_mu+slope_mu_Ui*Ui)*(longi_t+knots_start) + Ai*(inter_A+inter_A_Ui*Ui) - Ai*(slope_A+slope_A_Ui*Ui)*(longi_t+knots_start)) + sum(PsiX*Xi);
  NumericVector beta_1=(slope_mu + slope_mu_Ui*Ui + (slope_A+slope_A_Ui*Ui)*Ai);
  NumericVector alpha_vec = int2cpp((xb + xa)/2,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui,sty_end,surv_knots,surv_par_0,surv_par_1);
  
  double P=exp((alphaU*Ui)+alpha_vec[1]+(Ai*alpha_vec[2])+sum(Xi*alphaX));
  double mu_norm=(sum(beta_1*(yi-beta_0)) - tau*tau*P)/(sum(beta_1*beta_1));
  double sigma_norm=tau*tau/(sum(beta_1*beta_1));
  
  double cons_part= (0.5*mu_norm*mu_norm)/sigma_norm - sum((yi-beta_0)*(yi-beta_0))/(2*tau*tau) + alpha_vec[1] + alpha_vec[2]*Ai + int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui,sty_end,surv_knots,surv_par_0,surv_par_1)[0] - alpha_vec[3];
  
  // warning. Should apply for multiple surv_knots. 
  // if ((xb+xa)/2 > surv_knots[0]) {
  //   cons_part= (0.5*mu_norm*mu_norm)/sigma_norm - sum((yi-beta_0)*(yi-beta_0))/(2*tau*tau) + alpha_vec[1] + alpha_vec[2]*Ai + 
  //     (P-exp((alphaU*Ui)+alpha0+(Ai*alpha1)+sum(Xi*alphaX)))*surv_knots[0] + 
  //     int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui,sty_end,surv_knots,surv_par_0,surv_par_1)[0];
  // }
  
  double inte_result=sqrt(2*M_PI*sigma_norm) * (R::pnorm(xb,mu_norm,sqrt(sigma_norm),1,0)-R::pnorm(xa,mu_norm,sqrt(sigma_norm),1,0));
  // Rcout << " mu_norm " << mu_norm << " sigma_norm " << sqrt(sigma_norm) << " part1 " <<R::pnorm(xb,mu_norm,sigma_norm,1,0) << " part2 " << R::pnorm(xa,mu_norm,sigma_norm,1,0) << " const " << cons_part << "\n";
  // Rcout << " beta_1 " << beta_1 << " beta_0 " << beta_0 << " tstar_tend_comp " << tstar_tend_comp << "\n";
  // Rcout << " 0.5*mu_norm*mu_norm/sigma_norm " << 0.5*mu_norm*mu_norm << " sum((yi-beta_1)*(yi-beta_1))/(2*tau*tau) " << sum((yi-beta_1)*(yi-beta_1))/(2*tau*tau) << "\n";
  NumericVector result={inte_result,cons_part,mu_norm};
  return result;
}


NumericVector findmaxMtrmout3(double ti,NumericVector yi,double tau,double sigma,NumericVector PsiX,double gamma0,double gamma1,
                              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
                              double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
                              NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec,
                              NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, 
                              double pot_max, double mu_value, double pot_knots_start, double pot_knots_end, int left_right) {
  
  
  double left;
  double right;
  double M_temp;
  double M;
  double xa=0;
  double xb=0;
  double tol_length=1e-6;
  left=pot_knots_start;
  right=pot_knots_end;
  M=pot_max;
  // 0=left
  if (left_right==0) {
    xa=left;
    M_temp=trm3outcpp(right , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                      alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                      beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    if (M_temp>M-20) {
      xb=right;
    } else {
      M_temp=M+10;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>tol_length) {
        xb=(left+right)/2;
        M_temp=trm3outcpp(xb , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                          alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                          beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
        if (M_temp>M-20) {
          left=xb;
        } else {
          right=xb;
        }
      }
    }
  } else if (left_right==1){
    xb=right;
    M_temp=trm3outcpp(left , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                      alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                      beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    if (M_temp>M-20) {
      xa=left;
    } else {
      M_temp=M+10;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>tol_length) {
        xa=(left+right)/2;
        M_temp=trm3outcpp(xa , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                          alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                          beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
        if (M_temp>M-20) {
          right=xa;
        } else {
          left=xa;
        }
      }
    }
  } else {
    left=pot_knots_start;
    right=mu_value;
    
    M_temp=trm3outcpp(left , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                      alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                      beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    if (M_temp>M-20) {
      xa=left;
    } else {
      M_temp=M+10;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>tol_length) {
        xa=(left+right)/2;
        M_temp=trm3outcpp(xa , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                          alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                          beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
        if (M_temp>M-20) {
          right=xa;
        } else {
          left=xa;
        }
      }
    }
    
    
    left=mu_value;
    right=pot_knots_end;
    M_temp=trm3outcpp(right , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                      alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                      beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    if (M_temp>M-20) {
      xb=right;
    } else {
      M_temp=M+10;
      while(std::abs(M_temp-M+20)>0.001 && std::abs(left-right)>tol_length) {
        xb=(left+right)/2;
        M_temp=trm3outcpp(xb , ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                          alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                          beta_mu_vec, beta_A_vec,beta_mu_Ui_vec,beta_A_Ui_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
        if (M_temp>M-20) {
          left=xb;
        } else {
          right=xb;
        }
      }
    }
    // Rcout << " something wrong with inner find xa xb" << "\n";
  }
  NumericVector result={xa,xb};
  return(result);
}









double fftrm3cpp(double y, double xa, double xb, double ti,NumericVector yi,double tau,double sigma,NumericVector PsiX,double gamma0,double gamma1,
                 NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
                 double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
                 NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                 NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double M) {
  double u = -1; 
  double v = 1 ; 
  double duv = (xb-xa)/(v-u);
  // Rcout << " duv* " << exp(-M+trm3outcpp(xa + duv*(y-u), ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                                  alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 )) 
  //   << " duv " << duv << " duv*(y-u)" << xa+duv*(y-u)  << " expin " << -M+trm3outcpp(xa + duv*(y-u), ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                                  alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ) 
  //   << " M " << -M << " expin2 " << trm3outcpp(xa + duv*(y-u), ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                                  alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ) << "\n";
  return duv * exp(-M+trm3outcpp(xa + duv*(y-u), ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                 alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                                 knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ));
  
}

double quadinfIntcpptrm3(double xa, double xb, double tol, List nodes, List weights,double ti,NumericVector yi,double tau,double sigma,
                         NumericVector PsiX,double gamma0,double gamma1,
                         NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
                         double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
                         NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                         NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double M) {
  double Q;
  double h=0.5;
  NumericVector x=nodes[0];
  NumericVector w=weights[0];
  double s= w[6] * fftrm3cpp(x[6],xa, xb, ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                             alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t, beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                             knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M);
  for (int j=0;j<6;j++) {
    s=s + w[j] * (fftrm3cpp(x[j],xa, xb, ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t, beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M) + 
                              fftrm3cpp(-x[j],xa, xb, ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                        alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t, beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                                        knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M));
  }
  Q=s*h;
  //int k1;
  double delta1;
  double newQ;
  for (int k=1;k<7;k++) {
    x=nodes[k];
    w=weights[k];
    s=0;
    for (int j=0;j<w.size();j++) {
      s=s+w[j]*(fftrm3cpp(x[j],xa, xb, ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                          alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t, beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                          knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M) + 
                            fftrm3cpp(-x[j],xa, xb, ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                      alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t, beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec,
                                      knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M));
    }
    h = h/2;
    newQ=s*h + Q/2;
    delta1 = std::abs(newQ-Q);
    Q=newQ;
    if (delta1<tol) {
      break;
    }
  }
  
  //h=h/2;
  //newQ=s*h + Q/2;
  //NumericVector rest(4);
  //rest[0]=Q;
  //rest[1]=k1;
  //rest[2]=s;
  //rest[3]=h;
  
  return Q;
  
  //  if (xa==xb) {
  //    Q=0;
  //  } else if (xa>xb) {
  //    std::cout "xa should be smaller than xb"
  //  }
  
  
}

double li2cpp(double ti,NumericVector yi,double tau,NumericVector PsiX,double gamma0,double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
              double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
              NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec,
              NumericVector knots, double t_end, double sty_end, 
              NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double tol,
              List nodes,List weights){
  
  int ni=tstari.size();
  double sigma=1;
  double Ui_a=Ui;
  double trm1= -ni*log(tau);
  double trm21=Ui_a;
  double trm22=(alphaU+gammaU-Ui_a/(2*sigma*sigma));
  double trm2=trm21*trm22+sum(Xi*alphaX);
  // t_end should be the end of knots for integration
  // Rcout << " Ui " << Ui << " li2start " << "\n";
  NumericVector knots_2=clone(knots);
  if (std::abs(max(knots)-t_end)>1e-6) {
    knots_2.push_back(t_end);
  }
  
  
  int K=knots_2.size()-1;
  int L=longi_t.size();
  int a=0;
  // Rcout << " K " << K << " L " << L << " Ui " << Ui << "\n";
  NumericVector pot_knots_temp(L*(K+1));
  for (int i=0;i<L;i++) {
    for (int j=0;j<K+1;j++) {
      pot_knots_temp[a]=longi_t[i]+knots_2[j];
      a=a+1;
    }
  }
  
  for (int z=0;z<surv_knots.size();z++) {
    pot_knots_temp.push_back(surv_knots[z] + 1e-5);
    pot_knots_temp.push_back(surv_knots[z] - 1e-5);
  }
  
  // Rcout << " a " << a << "\n";
  NumericVector pot_knots_unsort=pot_knots_temp[pot_knots_temp>=ti];
  // Rcout << " Ui " << Ui << " li2start " << " pot_knots_temp " << pot_knots_temp << "\n";
  
  NumericVector pot_knots = sort_unique(pot_knots_unsort);
  // Rcout << " pot_knots " << pot_knots << " UI " << Ui << "\n";
  if (std::abs(pot_knots[0]-ti)>1e-6) {
    pot_knots.push_front(ti);
  }
  // Rcout << " pot_knots " << pot_knots << " UI " << Ui << "\n";
  
  int num_pot_knots=pot_knots.size();
  NumericVector pot_const(num_pot_knots-1);
  NumericVector pot_inte(num_pot_knots-1);
  NumericVector inte_temp;
  NumericVector pot_mu(num_pot_knots-1);
  LogicalVector not_na(num_pot_knots-1);
  NumericVector alpha_vec;
  for (int k=0;k<num_pot_knots-1;k++) {
    inte_temp=intergrate_trm3outcpp(pot_knots[k],pot_knots[k+1],ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                                    alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                                    beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1);
    
    pot_inte[k]=inte_temp[0];
    pot_const[k]=inte_temp[1];
    pot_mu[k]=inte_temp[2];
  }
  // Rcout << " Ui " << Ui << " li2start " << " pot_inte " << pot_inte << "\n";
  // Rcout << " pot_knots " << pot_knots << " pot_inte " << pot_inte << "\n";
  not_na=Rcpp::is_nan(pot_inte);
  // Rcout << " not_na " << not_na << "\n";
  pot_inte=pot_inte[!not_na];
  pot_const=pot_const[!not_na];
  pot_mu=pot_mu[!not_na];
  NumericVector pot_knots_ori=clone(pot_knots);
  NumericVector pot_knots_ori3=clone(pot_knots);
  pot_knots.erase(num_pot_knots-1);
  pot_knots_ori3.erase(0);
  // Rcout << " step1 "  << " Ui " << Ui << "\n";
  pot_knots=pot_knots[!not_na];
  pot_knots_ori3=pot_knots_ori3[!not_na];
  NumericVector pot_knots_ori2=clone(pot_knots);
  // pot_knots.push_back(pot_knots_ori[pot_knots.size()]);
  pot_knots.push_back(pot_knots_ori3[pot_knots.size()-1]);
  // NumericVector pot_const_ori=clone(pot_const);
  NumericVector pot_max(pot_knots.size()-1);
  IntegerVector left_right(pot_knots.size()-1);
  
  // Rcout << " Ui " << Ui << " li2start " << " pot_knots " << pot_knots << "\n";
  
  for (int f=0;f<pot_knots.size()-1;f++) {
    // inte_temp=intergrate_trm3outcpp(pot_knots[f],pot_knots[f+1],ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
    //                       alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
    //                       beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    // mu_temp[f] = pot_mu[f];
    if (pot_mu[f] > pot_knots[f+1]) {
      pot_max[f]=trm3outcpp(pot_knots[f+1],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=1;
    } else if (pot_mu[f] < pot_knots[f]) {
      pot_max[f]=trm3outcpp(pot_knots[f],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=0;
    } else {
      pot_max[f]=trm3outcpp(pot_mu[f],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=2;
    }
  }
  // Rcout << " Ui " << Ui << " li2start " << " pot_max " << pot_max << "\n";
  
  // int max_index=which_max(pot_max);
  // double max_inte_value=pot_inte[max_index];
  // double const_max=pot_const[max_index];
  double max_fun_value=max(pot_max);
  pot_inte=pot_inte[pot_max>=max_fun_value-20];
  pot_const=pot_const[pot_max>=max_fun_value-20];
  pot_knots_ori2=pot_knots_ori2[pot_max>=max_fun_value-20];
  pot_knots_ori3=pot_knots_ori3[pot_max>=max_fun_value-20];
  left_right=left_right[pot_max>=max_fun_value-20];
  NumericVector pot_max1=pot_max[pot_max>=max_fun_value-20];
  pot_mu=pot_mu[pot_max>=max_fun_value-20];
  
  
  LogicalVector no_too_close= (pot_knots_ori3 - pot_knots_ori2) > 1e-4;
  pot_inte=pot_inte[no_too_close];
  pot_const=pot_const[no_too_close];
  pot_knots_ori2=pot_knots_ori2[no_too_close];
  pot_knots_ori3=pot_knots_ori3[no_too_close];
  left_right=left_right[no_too_close];
  pot_max1=pot_max1[no_too_close];
  pot_mu=pot_mu[no_too_close];
  
  double trm3;
  // double trm31;
  double inner_last=trm3outcpp(pot_knots[pot_knots.size()-1],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                               alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                               beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                               knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
  alpha_vec = int2cpp(pot_knots[pot_knots.size()-1],alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1);
  double M;
  double xa;
  double xb;
  double knots_start_point;
  double knots_end_point;
  int left_right_status;
  NumericVector xa_xb_vec(2);
  // if (const_max<inner_last) {
  //   const_max=inner_last;
  // }
  // Rcout << " pot_max1 " << pot_max1 << " max_fun_value " << max_fun_value << "\n";
  // Rcout << " const_max " << const_max << " pot_const " << pot_const  << " pot_inte " << pot_inte << " pot_mu " << pot_mu << " pot_knots " << pot_knots << "\n";
  NumericVector fin_inte(pot_inte.size());
  NumericVector fin_max(pot_inte.size());
  
  double mu_value;
  
  // Rcout << " pot_inte " << pot_inte << " pot_const " << pot_const << " pot_knots_ori2 " << pot_knots_ori2 << " pot_knots_ori3 " << pot_knots_ori3 << " left_right " << left_right << " pot_max1 " << pot_max1 << " Ui "<< Ui << "\n";
  for (int p=0;p<pot_inte.size();p++) {
    if (pot_inte[p]>1e-3 || left_right[p]==2) {
      fin_inte[p]=pot_inte[p];
      fin_max[p]=pot_const[p];
    } else {
      M=pot_max1[p];
      knots_start_point=pot_knots_ori2[p];
      knots_end_point=pot_knots_ori3[p];
      left_right_status=left_right[p];
      mu_value=pot_mu[p];
      xa_xb_vec=findmaxMtrmout3(ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,
                                beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M ,mu_value, knots_start_point ,knots_end_point,left_right_status);
      xa=xa_xb_vec[0];
      xb=xa_xb_vec[1];
      // Rcout << " xa " << xa << " xb " << xb << " M " << M << "\n";
      // fin_inte[p]=quadinfIntcpptrm3(knots_start_point, knots_end_point, tol,nodes, weights,ti, yi, tau, sigma, 
      //                               PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
      //                               alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,
      //                               beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M);
      
      fin_inte[p]=quadinfIntcpptrm3(xa, xb, tol,nodes, weights,ti, yi, tau, sigma,
                                    PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                                    alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,
                                    beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M);
      
      
      fin_max[p]=M;
      // Rcout << " M " << M << " pot_knots_ori2[p] " << pot_knots_ori2[p] << " pot_knots_ori3[p] " << pot_knots_ori3[p] << " xa_xb_vec " << xa_xb_vec << " Ui " << Ui << "\n";
    }
  }
  // Rcout << " fin_inte " << fin_inte << " fin_max " << fin_max << " Ui "<< Ui << "\n";
  max_fun_value=max(fin_max);
  if (max_fun_value<inner_last) {
    max_fun_value = inner_last;
  }
  
  double inf_check=exp((inner_last-max_fun_value)-((alphaU*Ui)+alpha_vec[1]+(Ai*alpha_vec[2])+sum(Xi*alphaX)));
  trm3 = max_fun_value + log(sum(exp(fin_max-max_fun_value)*fin_inte) + inf_check);
  if (std::isinf(inf_check)) {
    trm3=max_fun_value + (inner_last-max_fun_value)-((alphaU*Ui)+alpha_vec[1]+(Ai*alpha_vec[2])+sum(Xi*alphaX));
  }
  // Rcout << " Ui " << Ui << " li2end " << "\n";
  
  
  
  
  
  // if (xa<xb1) {
  //   if (M<M1) {
  //     M=M1;
  //   }
  // 
  //   
  //   
  //   
  // } else {
  //   trm31=quadinfIntcpptrm3(xa, xb, tol,nodes, weights,ti, yi, tau, sigma, 
  //                           PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                           alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,M);
  //   trm32=quadinfIntcpptrm3(xa1, xb1, tol,nodes, weights,ti, yi, tau, sigma, 
  //                           PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                           alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,M1);
  //   
  //   M2=std::max(M,M1);
  //   trm3=M2+log(exp(M-M2)*trm31+exp(M1-M2)*trm32);
  // }
  
  //Rcout << "M " << M << " int " << quadinfIntcpptrm3(ti,tol,nodes, weights,ti, yi, tau, sigma, 
  //                                                 PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                                                 alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,M) << "\n" ;
  // Rcout << " Minner " << M << " xa " << xa << " xb " << xb << " trm31 " << trm31 << "\n";
  NumericVector gamma_vec=int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1);
  double trm4=gamma_vec[1];
  double trm5=gamma_vec[2]*Ai;
  double trm6=sum(Xi*gammaX);
  double trm7= -gamma_vec[0];
  double trm8= -log(sigma);
  double trm9 = -int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1)[0];
  double rtn2=trm1+trm2+trm3+trm4+trm5+trm6+trm7+trm8+trm9;
  return rtn2;
  
}


double li2cpp_shell(double Ui, List par_list) {
  double ti=par_list[0];
  NumericVector yi=par_list[1];
  NumericVector thetam=par_list[2];
  NumericVector PsiX=par_list[3];
  NumericVector Xi=par_list[4];
  double Ai=par_list[5];
  NumericVector tstari=par_list[6];
  NumericVector beta_mu_vec=par_list[7];
  NumericVector beta_A_vec=par_list[8];
  NumericVector knots=par_list[9];
  double t_end=par_list[10];
  double sty_end=par_list[11];
  NumericVector surv_knots=par_list[12];
  NumericVector drop_knots=par_list[13];
  NumericVector surv_par_0=par_list[14];
  NumericVector surv_par_1=par_list[15];
  NumericVector drop_par_0=par_list[16];
  NumericVector drop_par_1=par_list[17];
  
  
  List add_elements=par_list[18];
  
  
  
  NumericVector longi_t=add_elements[0];
  double tol=add_elements[1];
  List nodes=add_elements[2];
  List weights=add_elements[3];
  
  NumericVector alphaX=add_elements[4];
  NumericVector gammaX =add_elements[5];
  NumericVector beta_mu_Ui_vec=add_elements[6];
  NumericVector beta_A_Ui_vec=add_elements[7];
  
  double betamu=thetam[0];
  double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];
  
  double alphaU=thetam[4];
  double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];
  
  double gammaU=thetam[8];

  double res=li2cpp(ti,yi,tau,PsiX,gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, betamu,
                    betaA, Xi, Ai, tstari, Ui,longi_t,beta_mu_vec, beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, 
                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ,tol,nodes,weights);
  return res;
}




double li3cpp(double ti,NumericVector yi,double tau,NumericVector PsiX,double gamma0,double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,double betamu,
              double betaA,NumericVector Xi,double Ai,NumericVector tstari,double Ui,NumericVector longi_t,
              NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
              NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1,double tol,
              List nodes,List weights){
  
  double sigma=1;
  int ni=tstari.size();
  double Ui_a=Ui;
  // Rcout << " Ui " << Ui << "\n";
  double trm1= -ni*log(tau);
  double trm2 = ((alphaU-Ui_a/(2*sigma*sigma))*Ui_a)+sum(Xi*alphaX);
  // Rcout << " step 3 " << "\n";
  // NumericVector Ms={-2416.87,31.8722,34.4985};
  
  
  
  // Integration
  NumericVector knots_2=clone(knots);
  if (std::abs(max(knots)-t_end)>1e-6) {
    knots_2.push_back(t_end);
  }
  
  
  int K=knots_2.size()-1;
  int L=longi_t.size();
  int a=0;
  // Rcout << " K " << K << " L " << L << " Ui " << Ui << "\n";
  NumericVector pot_knots_temp(L*(K+1));
  for (int i=0;i<L;i++) {
    for (int j=0;j<K+1;j++) {
      pot_knots_temp[a]=longi_t[i]+knots_2[j];
      a=a+1;
    }
  }
  
  for (int z=0;z<surv_knots.size();z++) {
    pot_knots_temp.push_back(surv_knots[z] + 1e-6);
    pot_knots_temp.push_back(surv_knots[z] - 1e-6);
  }
  
  // Rcout << " a " << a << "\n";
  NumericVector pot_knots_unsort=pot_knots_temp[pot_knots_temp>=ti];
  // Rcout << " Ui " << Ui << " li2start " << " pot_knots_temp " << pot_knots_temp << "\n";
  
  NumericVector pot_knots = sort_unique(pot_knots_unsort);
  // Rcout << " pot_knots " << pot_knots << " UI " << Ui << "\n";
  if (std::abs(pot_knots[0]-ti)>1e-6) {
    pot_knots.push_front(ti);
  }
  // Rcout << " pot_knots " << pot_knots << " UI " << Ui << "\n";
  
  int num_pot_knots=pot_knots.size();
  NumericVector pot_const(num_pot_knots-1);
  NumericVector pot_inte(num_pot_knots-1);
  NumericVector inte_temp;
  NumericVector pot_mu(num_pot_knots-1);
  LogicalVector not_na(num_pot_knots-1);
  NumericVector alpha_vec;
  for (int k=0;k<num_pot_knots-1;k++) {
    inte_temp=intergrate_trm3outcpp(pot_knots[k],pot_knots[k+1],ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                                    alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                                    beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1);
    
    pot_inte[k]=inte_temp[0];
    pot_const[k]=inte_temp[1];
    pot_mu[k]=inte_temp[2];
  }
  // Rcout << " Ui " << Ui << " li2start " << " pot_inte " << pot_inte << "\n";
  // Rcout << " pot_knots " << pot_knots << " pot_inte " << pot_inte << "\n";
  not_na=Rcpp::is_nan(pot_inte);
  // Rcout << " not_na " << not_na << "\n";
  pot_inte=pot_inte[!not_na];
  pot_const=pot_const[!not_na];
  pot_mu=pot_mu[!not_na];
  NumericVector pot_knots_ori=clone(pot_knots);
  NumericVector pot_knots_ori3=clone(pot_knots);
  pot_knots.erase(num_pot_knots-1);
  pot_knots_ori3.erase(0);
  // Rcout << " step1 "  << " Ui " << Ui << "\n";
  pot_knots=pot_knots[!not_na];
  pot_knots_ori3=pot_knots_ori3[!not_na];
  NumericVector pot_knots_ori2=clone(pot_knots);
  // pot_knots.push_back(pot_knots_ori[pot_knots.size()]);
  pot_knots.push_back(pot_knots_ori3[pot_knots.size()-1]);
  
  // NumericVector pot_const_ori=clone(pot_const);
  NumericVector pot_max(pot_knots.size()-1);
  IntegerVector left_right(pot_knots.size()-1);
  
  // Rcout << " Ui " << Ui << " li2start " << " pot_knots " << pot_knots << "\n";
  
  for (int f=0;f<pot_knots.size()-1;f++) {
    // inte_temp=intergrate_trm3outcpp(pot_knots[f],pot_knots[f+1],ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
    //                       alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
    //                       beta_mu_vec, beta_A_vec, knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
    // mu_temp[f] = pot_mu[f];
    if (pot_mu[f] > pot_knots[f+1]) {
      pot_max[f]=trm3outcpp(pot_knots[f+1],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=1;
    } else if (pot_mu[f] < pot_knots[f]) {
      pot_max[f]=trm3outcpp(pot_knots[f],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=0;
    } else {
      pot_max[f]=trm3outcpp(pot_mu[f],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                            alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                            beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                            knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
      left_right[f]=2;
    }
  }
  // Rcout << " Ui " << Ui << " li2start " << " pot_max " << pot_max << "\n";
  
  // int max_index=which_max(pot_max);
  // double max_inte_value=pot_inte[max_index];
  // double const_max=pot_const[max_index];
  double max_fun_value=max(pot_max);
  pot_inte=pot_inte[pot_max>=max_fun_value-20];
  pot_const=pot_const[pot_max>=max_fun_value-20];
  pot_knots_ori2=pot_knots_ori2[pot_max>=max_fun_value-20];
  pot_knots_ori3=pot_knots_ori3[pot_max>=max_fun_value-20];
  left_right=left_right[pot_max>=max_fun_value-20];
  NumericVector pot_max1=pot_max[pot_max>=max_fun_value-20];
  pot_mu=pot_mu[pot_max>=max_fun_value-20];
  
  
  LogicalVector no_too_close= (pot_knots_ori3 - pot_knots_ori2) > 1e-4;
  pot_inte=pot_inte[no_too_close];
  pot_const=pot_const[no_too_close];
  pot_knots_ori2=pot_knots_ori2[no_too_close];
  pot_knots_ori3=pot_knots_ori3[no_too_close];
  left_right=left_right[no_too_close];
  pot_max1=pot_max1[no_too_close];
  pot_mu=pot_mu[no_too_close];
  
  double trm3;
  // double trm31;
  double inner_last=trm3outcpp(pot_knots[pot_knots.size()-1],ti, yi, tau, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1,
                               alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,
                               beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                               knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 );
  alpha_vec = int2cpp(pot_knots[pot_knots.size()-1],alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1);
  double M;
  double xa;
  double xb;
  double knots_start_point;
  double knots_end_point;
  int left_right_status;
  NumericVector xa_xb_vec(2);
  // if (const_max<inner_last) {
  //   const_max=inner_last;
  // }
  // Rcout << " pot_max1 " << pot_max1 << " max_fun_value " << max_fun_value << "\n";
  // Rcout << " const_max " << const_max << " pot_const " << pot_const  << " pot_inte " << pot_inte << " pot_mu " << pot_mu << " pot_knots " << pot_knots << "\n";
  NumericVector fin_inte(pot_inte.size());
  NumericVector fin_max(pot_inte.size());
  
  double mu_value;
  
  // Rcout << " pot_inte " << pot_inte << " pot_const " << pot_const << " pot_knots_ori2 " << pot_knots_ori2 << " pot_knots_ori3 " << pot_knots_ori3 << " left_right " << left_right << " pot_max1 " << pot_max1 << " Ui "<< Ui << "\n";
  for (int p=0;p<pot_inte.size();p++) {
    if (pot_inte[p]>1e-3 || left_right[p]==2) {
      fin_inte[p]=pot_inte[p];
      fin_max[p]=pot_const[p];
    } else {
      M=pot_max1[p];
      knots_start_point=pot_knots_ori2[p];
      knots_end_point=pot_knots_ori3[p];
      left_right_status=left_right[p];
      mu_value=pot_mu[p];
      
      xa_xb_vec=findmaxMtrmout3(ti, yi, tau, sigma, PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,
                                beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M ,mu_value, knots_start_point ,knots_end_point,left_right_status);
      xa=xa_xb_vec[0];
      xb=xa_xb_vec[1];
      // Rcout << " xa " << xa << " xb " << xb << " M " << M << "\n";
      fin_inte[p]=quadinfIntcpptrm3(xa, xb, tol,nodes, weights,ti, yi, tau, sigma, 
                                    PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
                                    alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui_a, longi_t,
                                    beta_mu_vec, beta_A_vec, beta_mu_Ui_vec,beta_A_Ui_vec, 
                                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 , M);
      fin_max[p]=M;
      // Rcout << " M " << M << " pot_knots_ori2[p] " << pot_knots_ori2[p] << " pot_knots_ori3[p] " << pot_knots_ori3[p] << " xa_xb_vec " << xa_xb_vec << " Ui " << Ui << "\n";
    }
  }
  // Rcout << " fin_inte " << fin_inte << " fin_max " << fin_max << " Ui "<< Ui << "\n";
  max_fun_value=max(fin_max);
  if (max_fun_value<inner_last) {
    max_fun_value = inner_last;
  }
  double inf_check=exp((inner_last-max_fun_value)-((alphaU*Ui)+alpha_vec[1]+(Ai*alpha_vec[2])+sum(Xi*alphaX)));
  trm3 = max_fun_value + log(sum(exp(fin_max-max_fun_value)*fin_inte) + inf_check);
  if (std::isinf(inf_check)) {
    trm3=max_fun_value + (inner_last-max_fun_value)-((alphaU*Ui)+alpha_vec[1]+(Ai*alpha_vec[2])+sum(Xi*alphaX));
  }
  
  // Rcout << " Ui " << Ui << " li2end " << "\n";
  
  // Rcout << " fin_max " << fin_max << " pot_inte " << pot_inte << " fin_inte " << fin_inte << " max_fun_value " << max_fun_value << " pot_const " << pot_const << "\n";
  
  double trm4= - int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1)[0];
  double trm5= - log(sigma);
  double trm6 = -int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1)[0];
  double rtn = trm1+trm2+trm3+trm4+trm5+trm6;
  // Rcout << " rtn " << rtn << " Ui " << Ui << "\n";
  // Rcout << " trm1 " << trm1 << " trm2 " << trm2 << " trm3 " << trm3 << " trm4 " << trm4 << " trm5 " << trm5 << " trm6 " << trm6 << "\n";
  // Rcout << " trm1 " << trm1  << " trm3 " << trm3 << " trm4 " << trm4  << " trm6 " << trm6 << "\n";
  // Rcout << " Ui " << Ui << " trm3 " << quadinfIntcpptrm3(ti,tol,nodes, weights,ti, yi, tau, sigma, 
  //                                                    PsiX, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, 
  //                                                    alphaX, alphaU, betamu, betaA, Xi, Ai, tstari, Ui, longi_t,M) << " innerM " << M << "\n";
  // Rcout << " Ui " << Ui << " trm31 " << trm31 << " trm3 " << trm3 << " rtn " <<rtn <<  "\n";
  return rtn;
  
  
}

double li3cpp_shell(double Ui, List par_list) {
  double ti=par_list[0];
  NumericVector yi=par_list[1];
  NumericVector thetam=par_list[2];
  NumericVector PsiX=par_list[3];
  NumericVector Xi=par_list[4];
  double Ai=par_list[5];
  NumericVector tstari=par_list[6];
  NumericVector beta_mu_vec=par_list[7];
  NumericVector beta_A_vec=par_list[8];
  NumericVector knots=par_list[9];
  double t_end=par_list[10];
  double sty_end=par_list[11];
  NumericVector surv_knots=par_list[12];
  NumericVector drop_knots=par_list[13];
  NumericVector surv_par_0=par_list[14];
  NumericVector surv_par_1=par_list[15];
  NumericVector drop_par_0=par_list[16];
  NumericVector drop_par_1=par_list[17];
  
  List add_elements=par_list[18];
  

  
  NumericVector longi_t=add_elements[0];
  double tol=add_elements[1];
  List nodes=add_elements[2];
  List weights=add_elements[3];
  NumericVector alphaX=add_elements[4];
  NumericVector gammaX =add_elements[5];
  NumericVector beta_mu_Ui_vec=add_elements[6];
  NumericVector beta_A_Ui_vec=add_elements[7];
  
  double betamu=thetam[0];
  double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];
  
  double alphaU=thetam[4];
  double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];
  
  double gammaU=thetam[8];

  double res=li3cpp(ti,yi,tau,PsiX,gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, betamu,
                    betaA, Xi, Ai, tstari, Ui,longi_t,beta_mu_vec, beta_A_vec, beta_mu_Ui_vec, beta_A_Ui_vec, 
                    knots, t_end, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1 ,tol,nodes,weights);
  return res;
}

double li4cpp(double ti, double gamma0, double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,
              NumericVector Xi,double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double Ui) {
  double sigma=1;
  double Ui_a=Ui;
  double trm1 = -int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1)[0];
  NumericVector alpha_vec=int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1);
  double trm2=alpha_vec[1];
  double trm3=alpha_vec[2]*Ai;
  double trm4=(alphaU-Ui_a/(2*sigma*sigma))*Ui_a;
  double trm5=sum(Xi*alphaX);
  double trm6= -alpha_vec[0];
  double trm7= -log(sigma);
  double ret=trm1+trm2+trm3+trm4+trm5+trm6+trm7;
  return ret;
}

double li4cpp_shell(double Ui, List par_list) {
  
  double ti=par_list[0];
  // NumericVector yi=par_list[1];
  NumericVector thetam=par_list[1];
  NumericVector PsiX=par_list[2];
  NumericVector Xi=par_list[3];
  double Ai=par_list[4];
  // NumericVector tstari=par_list[6];
  // NumericVector beta_mu_vec=par_list[7];
  // NumericVector beta_A_vec=par_list[8];
  // NumericVector knots=par_list[9];
  // double t_end=par_list[10];
  double sty_end=par_list[5];
  NumericVector surv_knots=par_list[6];
  NumericVector drop_knots=par_list[7];
  NumericVector surv_par_0=par_list[8];
  NumericVector surv_par_1=par_list[9];
  NumericVector drop_par_0=par_list[10];
  NumericVector drop_par_1=par_list[11];
  
  NumericVector alphaX=par_list[12];
  NumericVector gammaX =par_list[13];
  
  // double betamu=thetam[0];
  // double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];
  
  double alphaU=thetam[4];
  // double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];
  
  double gammaU=thetam[8];

  double res=li4cpp(ti, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, 
                    Xi, Ai, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1, Ui);
  return res;
}

double li5cpp(double ti, double gamma0, double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,
              NumericVector Xi,double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double Ui){
  double sigma=1;
  double Ui_a=Ui;
  NumericVector gamma_vec=int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1);
  double trm1= -int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1)[0];
  double trm2=gamma_vec[1];
  double trm3=Ai*gamma_vec[2];
  double trm4=(gammaU-Ui_a/(2*sigma*sigma))*Ui_a;
  double trm5=sum(Xi*gammaX);
  double trm6= -gamma_vec[0];
  double trm7= -log(sigma);
  double rtn=trm1+trm2+trm3+trm4+trm5+trm6+trm7;
  return (rtn);
}

double li5cpp_shell(double Ui, List par_list) {
  double ti=par_list[0];
  // NumericVector yi=par_list[1];
  NumericVector thetam=par_list[1];
  NumericVector PsiX=par_list[2];
  NumericVector Xi=par_list[3];
  double Ai=par_list[4];
  // NumericVector tstari=par_list[6];
  // NumericVector beta_mu_vec=par_list[7];
  // NumericVector beta_A_vec=par_list[8];
  // NumericVector knots=par_list[9];
  // double t_end=par_list[10];
  double sty_end=par_list[5];
  NumericVector surv_knots=par_list[6];
  NumericVector drop_knots=par_list[7];
  NumericVector surv_par_0=par_list[8];
  NumericVector surv_par_1=par_list[9];
  NumericVector drop_par_0=par_list[10];
  NumericVector drop_par_1=par_list[11];
  
  NumericVector alphaX=par_list[12];
  NumericVector gammaX =par_list[13];
  
  // double betamu=thetam[0];
  // double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];
  
  double alphaU=thetam[4];
  // double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];
  
  double gammaU=thetam[8];

  double res=li5cpp(ti, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, 
                    Xi, Ai, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1, Ui);
  return res;
}

double li6cpp(double ti, double gamma0, double gamma1,
              NumericVector gammaX,double gammaU,double alpha0,double alpha1,NumericVector alphaX,double alphaU,
              NumericVector Xi,double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double Ui){
  double sigma=1;
  double Ui_a=Ui;
  double trm1= -int2cpp(ti,alpha0,alpha1,alphaX,alphaU,Xi,Ai,Ui_a,sty_end,surv_knots,surv_par_0,surv_par_1)[0];
  double trm2= -int1cpp(ti,gamma0,gamma1,gammaX,gammaU,Xi,Ai,Ui_a,sty_end,drop_knots,drop_par_0,drop_par_1)[0];
  double trm3= -log(sigma);
  double trm4= -Ui_a*Ui_a/(2*sigma*sigma);
  double rtn=trm1+trm2+trm3+trm4;
  return (rtn);
}

double li6cpp_shell(double Ui, List par_list) {
  double ti=par_list[0];
  // NumericVector yi=par_list[1];
  NumericVector thetam=par_list[1];
  NumericVector PsiX=par_list[2];
  NumericVector Xi=par_list[3];
  double Ai=par_list[4];
  // NumericVector tstari=par_list[6];
  // NumericVector beta_mu_vec=par_list[7];
  // NumericVector beta_A_vec=par_list[8];
  // NumericVector knots=par_list[9];
  // double t_end=par_list[10];
  double sty_end=par_list[5];
  NumericVector surv_knots=par_list[6];
  NumericVector drop_knots=par_list[7];
  NumericVector surv_par_0=par_list[8];
  NumericVector surv_par_1=par_list[9];
  NumericVector drop_par_0=par_list[10];
  NumericVector drop_par_1=par_list[11];
  
  NumericVector alphaX=par_list[12];
  NumericVector gammaX =par_list[13];
  
  
  // double betamu=thetam[0];
  // double betaA=thetam[1];
  
  double alpha0=thetam[2];
  double alpha1=thetam[3];
  
  double alphaU=thetam[4];
  // double tau=thetam[5];
  double gamma0=thetam[6];
  double gamma1=thetam[7];
  
  double gammaU=thetam[8];

  double res=li6cpp(ti, gamma0, gamma1, gammaX, gammaU, alpha0, alpha1, alphaX, alphaU, 
                    Xi, Ai, sty_end, surv_knots, drop_knots, surv_par_0, surv_par_1, drop_par_0, drop_par_1, Ui);
  return res;
}



// [[Rcpp::export]]
double Egudenintg1hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector yi, NumericVector Xi, double Ai, NumericVector tstari, 
                          NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec,
                          NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, 
                          NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double tol, List nodes,List weights) {
  List add_elements = Rcpp::List::create(alphaX,gammaX,beta_mu_Ui_vec,beta_A_Ui_vec);
  List par_list = Rcpp::List::create(ti,yi,thetam,PsiX,Xi,Ai,tstari,beta_mu_vec,beta_A_vec,knots,t_end,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,add_elements);
  
  NumericMatrix Ms=find_max(li1cpp_shell,par_list,-10,10);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li1cpp_shell, par_list, M[i], xa[i], xb[i]);
  }
  double M_max=max(M);
  // Rcout << " M_max " << M_max << "\n";
  
  double result;
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
}


// [[Rcpp::export]]
double Egudenintg2hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector yi, NumericVector Xi, double Ai, NumericVector tstari,NumericVector longi_t, 
                          NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                          NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double tol, List nodes,List weights) {

  List add_list = Rcpp::List::create(longi_t,tol,nodes,weights,alphaX,gammaX,beta_mu_Ui_vec,beta_A_Ui_vec);
  List par_list = Rcpp::List::create(ti,yi,thetam,PsiX,Xi,Ai,tstari,beta_mu_vec,beta_A_vec,knots,t_end,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,add_list);
  
  NumericMatrix Ms=find_max(li2cpp_shell,par_list,-3,3);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  // double denting3=0;
  
  
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li2cpp_shell, par_list, M[i], xa[i], xb[i]);
    while (std::isinf(inte_res[i])) {
      M[i]= M[i]+400;
      inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li2cpp_shell, par_list, M[i], xa[i], xb[i]);
      //Rcout << " denting3 " << denting3 << "\n";
      Rcout << "Something wrong with M in li2" << "\n";
      
      
    }
  }
  double M_max=max(M);
  double result;
  // Rcout << " M_max " << M_max << "\n";
  
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
}


// [[Rcpp::export]]
double Egudenintg3hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector yi, NumericVector Xi, double Ai, NumericVector tstari,NumericVector longi_t, 
                          NumericVector beta_mu_vec, NumericVector beta_A_vec, NumericVector beta_mu_Ui_vec, NumericVector beta_A_Ui_vec, 
                          NumericVector knots, double t_end, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, double tol, List nodes,List weights) {

  List add_list = Rcpp::List::create(longi_t,tol,nodes,weights,alphaX,gammaX,beta_mu_Ui_vec,beta_A_Ui_vec);
  List par_list = Rcpp::List::create(ti,yi,thetam,PsiX,Xi,Ai,tstari,beta_mu_vec,beta_A_vec,knots,t_end,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,add_list);

  NumericMatrix Ms=find_max(li3cpp_shell,par_list,-3,3);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  // double denting3=0;
  
  
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li3cpp_shell, par_list, M[i], xa[i], xb[i]);
    
    while (std::isinf(inte_res[i])) {
      M[i]= M[i]+400;
      inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li3cpp_shell, par_list, M[i], xa[i], xb[i]);
      
      //Rcout << " denting3 " << denting3 << "\n";
      Rcout << "Something wrong with M in li3" << "\n";
      
      
    }
  }
  double M_max=max(M);
  double result;
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
}

// [[Rcpp::export]]
double Egudenintg4hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector Xi, double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, 
                          double tol, List nodes,List weights) {
  
  List par_list = Rcpp::List::create(ti,thetam,PsiX,Xi,Ai,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,alphaX,gammaX);
  
  NumericMatrix Ms=find_max(li4cpp_shell,par_list,-10,10);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li4cpp_shell, par_list, M[i], xa[i], xb[i]);
  }
  double M_max=max(M);
  double result;
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
  
}


// [[Rcpp::export]]
double Egudenintg5hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector Xi, double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, 
                          double tol, List nodes,List weights) {

  List par_list = Rcpp::List::create(ti,thetam,PsiX,Xi,Ai,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,alphaX,gammaX);
  
  NumericMatrix Ms=find_max(li5cpp_shell,par_list,-10,10);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li5cpp_shell, par_list, M[i], xa[i], xb[i]);
  }
  double M_max=max(M);
  double result;
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
  
}

// [[Rcpp::export]]
double Egudenintg6hesscpp(NumericVector thetam, NumericVector alphaX, NumericVector gammaX, NumericVector PsiX, double ti, NumericVector Xi, double Ai, double sty_end, NumericVector surv_knots, NumericVector drop_knots, NumericVector surv_par_0, NumericVector surv_par_1, NumericVector drop_par_0, NumericVector drop_par_1, 
                          double tol, List nodes,List weights) {
  
  
  List par_list = Rcpp::List::create(ti,thetam,PsiX,Xi,Ai,sty_end,surv_knots,
                                     drop_knots,surv_par_0,surv_par_1,drop_par_0,drop_par_1,alphaX,gammaX);
  
  NumericMatrix Ms=find_max(li6cpp_shell,par_list,-10,10);
  
  NumericVector M=Ms (_,0);
  NumericVector xa=Ms(_,1);
  NumericVector xb=Ms(_,2);
  NumericVector inte_res(M.size());
  for (int i=0;i<M.size();i++) {
    inte_res[i]=quadinfIntcppdenfunc(tol,nodes, weights, li6cpp_shell, par_list, M[i], xa[i], xb[i]);
  }
  double M_max=max(M);
  double result;
  result = M_max + log(sum(exp(M-M_max)*inte_res));
  return (result);
  
}



