### Figure 2, TB mortality reductions 
# running out 200 years, probably unnecessarily long
## Rapid Molecular Test 
# 1. standard of care sensitivity (smear neg) 
std_smneg <- c(tbdx(sn_smneg=0.50)$mort_tb[10001],tbdx(sn_smneg=0.05)$mort_tb[10001])
rmt_base <- tbdx(sn_smpos=0.98,sn_smneg=0.50)$mort_tb[10001]
round(100*(1- rmt_base/std_smneg), 1)

# 2. intervention sensitivity (smear neg)
std_base <- tbdx()$mort_tb[10001]
rmt_smneg <- c(tbdx(sn_smpos=0.98,sn_smneg=0.25)$mort_tb[10001]
			,tbdx(sn_smpos=0.98,sn_smneg=0.75)$mort_tb[10001])
round(100*(1- rmt_smneg/std_base), 1)

# 3. standard of care sensitivity (smear pos) 
std_smpos <- c(tbdx(sn_smpos=1.0)$mort_tb[10001],tbdx(sn_smpos=0.50)$mort_tb[10001])
rmt_base <- tbdx(sn_smpos=0.98,sn_smneg=0.50)$mort_tb[10001]
round(100*(1- rmt_base/std_smpos), 1)

# 4. intervention sensitivity (smear pos)
std_base <- tbdx()$mort_tb[10001]
rmt_smpos <- c(tbdx(sn_smpos=0.80,sn_smneg=0.50)$mort_tb[10001],tbdx(sn_smpos=1.0,sn_smneg=0.50)$mort_tb[10001])
round(100*(1- rmt_smpos/std_base), 1)

## TB Culture
# 5. standard of care sensitivity (smear neg) 
std_smneg <- c(tbdx(sn_smneg=0.50)$mort_tb[10001],tbdx(sn_smneg=0.05)$mort_tb[10001])
tbc_base <- tbdx(dx_neg=3,lfu=0.2,sn_smpos=0.98,sn_smneg=0.85)$mort_tb[10001]
round(100*(1- tbc_base/std_smneg), 1)

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
