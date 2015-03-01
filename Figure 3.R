### Figure 3, TB mortality reductions 
# running out 200 years, probably unnecessarily long
## Rapid Molecular Test 
# 1. Proportion of TB patients accessing clinics 
std_base_1 <- tbdx()$mort_tb[10001]
rmt_smneg_1 <- c(tbdx(sn_smpos=0.98,sn_smneg=0.5,acc_neg = 0.53, acc_hiv=0.53, acc_aids = 0.53)$mort_tb[10001]
			,tbdx(sn_smpos=0.98,sn_smneg=0.5,acc_neg = 0.734,acc_hiv=0.734, acc_aids = 0.734)$mort_tb[10001])
round(100*(1- rmt_smneg_1/std_base_1), 1)


# 2. Relative infectiousness of smear-negative TB 
std_base_2 <- tbdx()$mort_tb[10001]
rmt_smneg_2 <- c(tbdx(sn_smpos=0.98,sn_smneg=0.5,i_smneg=0.42)$mort_tb[10001]
			,tbdx(sn_smpos=0.98,sn_smneg=0.5,i_smneg=0.02)$mort_tb[10001])
round(100*(1- rmt_smneg_2/std_base_2), 1)



##Tuberculosis culture 

# 3. Proportion of TB patients accessing clinics 
std_base_3 <- tbdx()$mort_tb[10001]
rmt_smneg_3 <- c(tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85, acc_neg = 0.53, acc_hiv=0.53, acc_aids = 0.53)$mort_tb[10001]
			,tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85,acc_neg = 0.93,acc_hiv=0.93, acc_aids = 0.93)$mort_tb[10001])
round(100*(1- rmt_smneg_3/std_base_3), 1)


# 4. Relative infectiousness of smear-negative TB 
std_base_4 <- tbdx()$mort_tb[10001]
rmt_smneg_4 <- c(tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85,i_smneg=0.42)$mort_tb[10001]
			,tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85,i_smneg=0.02)$mort_tb[10001])
round(100*(1- rmt_smneg_4/std_base_4), 1)




## Community ACF
std_base_5 <- tbdx()$mort_tb[10001]
rmt_smneg_5 <- c(tbdx(acc_neg = 0.53, acc_hiv=0.53, acc_aids = 0.53)$mort_tb[10001]
			,tbdx(acc_neg = 0.93,acc_hiv=0.93, acc_aids = 0.93)$mort_tb[10001])
round(100*(1- rmt_smneg_5/std_base_5), 1)






# 6. intervention sensitivity (smear pos)
std_smpos <- c(tbdx(sn_smpos=1.0)$mort_tb[10001],tbdx(sn_smpos=0.50)$mort_tb[10001])
tbc_base <- tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001]
round(100*(1- tbc_base/std_smpos), 1)

# 7. proportion of patients returning for results ( = 1-lfu )
std_base <- tbdx()$mort_tb[10001]
tbc_lfu <- c(tbdx(dx_neg=3,lfu=1-0.5,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001]
			, tbdx(dx_neg=3,lfu=1-1,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001])
round(100*(1- tbc_lfu/std_base), 1)

# 8. intervention sensitivity (smear neg)
std_base <- tbdx()$mort_tb[10001]
tbc_smneg <- c(tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.50)$mort_tb[10001]
			, tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=1.0)$mort_tb[10001]) 
round(100*(1- tbc_smneg/std_base), 1)

# 9. diagnostic delay for results (days), assume 30 days per month 
std_base <- tbdx()$mort_tb[10001]
tbc_dx <- c(tbdx(dx_neg=1/(1/4 + 60/360),lfu=0.2,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001]
			, tbdx(dx_neg=1/(1/4 + 10/360),lfu=0.2,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001]) 
round(100*(1- tbc_dx/std_base), 1)

# 10. intervention sensitivity (smear pos)
std_base <- tbdx()$mort_tb[10001]
tbc_smpos <- c(tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.8,sn_smneg=0.85)$mort_tb[10001]
			, tbdx(dx_neg=3,lfu=0.2,sn_smpos=1.0,sn_smneg=0.85)$mort_tb[10001]) 
round(100*(1- tbc_smpos/std_base), 1)
