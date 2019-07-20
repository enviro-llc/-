# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 22:45:58 2019

@author: LiShu
"""
'''
# HJ 25.3-2014 附录F污染物扩散迁移推荐模型
class F:
    def F2_4(ρ_b,ρ_s,P_ws): #见G.1
        θ = 1-ρ_b/ρ_s
        ρ_w = 1 # unit kg/L
        θ_ws = ρ_b*P_ws/ρ_w
        θ_as = θ-θ_ws
        return(θ,θ_ws,θ_as)
        
    def F1(θ,θ_ws,θ_as,D_a):
        Deff_s = D_a*θ_as**3.33/θ**2+D_w*θ_ws**3.33/(H*θ**2)
        
a = F
'''

# HJ 25.3-2014 附录A 暴露评估推荐模型

# 暴露途径1 经口摄入土壤途径
def OISER(OSIR_c,ED_c,EF_c,BW_c,OSIR_a,ED_a,EF_a,BW_a,ABS_0,AT_ca,AT_nc):
    #敏感用地致癌效应
    OISER_ca = 10**(-6)*(((OSIR_c*ED_c*EF_c)/BW_c)+((OSIR_a*ED_a*EF_a)/BW_a)*ABS_0)/AT_ca 
    #敏感用地非致癌效应
    OISER_nc = 10**(-6)*((OSIR_c*ED_c*EF_c*ABS_0)/(BW_c*AT_nc))
    return(OISER_ca,OISER_nc)
    
    
# 暴露途径2 皮肤接触土壤
def DCSER():
    # 儿童与成人暴露皮肤表面积
    SAE_c = 239*(H_c**0.417)*(BW_c**0.517)*SER_c
    SAE_a = 239*(H_a**0.417)*(BW_a**0.517)*SER_a
     #敏感用地致癌效应
    DCSER_ca = 10**(-6)*((SAE_c*SSAR_c*EF_c*ED_c*E_v*ABS_d)/(BW_c*AT_ca)+(SAE_a*SSAR_a*EF_a*ED_a*E_v*ABS_d)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    DCSER_nc = 10**(-6)*(SAE_c*SSAR_c*EF_c*ED_c*E_v*ABS_d)/(BW_c*AT_nc)
    return(DCSER_ca,DCSER_nc)


# 暴露途径3 吸入土壤颗粒物
def PISER():
    #敏感用地致癌效应(怀疑文件中原公式（A.7）有错？第一项的被除数应为BW_c？)
    PISER_ca = 10**(-6)*((PM_10*DAIR_c*ED_c*PIAF*(fspo*EFO_c+fspi*EFI_c)/(BW_a*AT_ca))+(PM_10*DAIR_a*ED_a*PIAF*(fspo*EFO_a+fspi*EFI_a)/(BW_a*AT_ca)))
    #敏感用地非致癌效应
    PISER_nc = 10**(-6)*(PM_10*DAIR_c*ED_c*PIAF*(fspo*EFO_c+fspi*EFI_c)/(BW_c*AT_nc))
    return(PISER_ca,PISER_nc)
    
    
# 暴露途径4 吸入室外空气中来自表层土壤的气态污染物
def IOVER1():
    #非饱和土层土壤中总孔隙体积比(F.2)、孔隙水体积比(F.3)、孔隙空气体积比(F.4)
    θ = 1-(ρ_b/ρ_s)   #ρ见G.1
    θ_ws = ρ_b*ρ_ws/ρ_w
    θ_as = θ-θ_ws 
    # 土壤中气态污染物的有效扩散系数(F.1)
    Deff_s = D_a*(θ_as**3.33)/(θ**2)+D_w*(θ_ws**3.33)/(H*θ**2) # D,H见B.2
    
    #土壤有机碳/土壤孔隙水分配系数
    K_d = K_oc*f_om/1700 #f_om为土壤有机质含量，需根据场地调查获得
    # 土壤-水污染物分配系数
    K_sw = (θ_ws+(K_d*ρ_b)+(H*θ_as))/ρ_b
    
    #室外空气中气态污染物扩散因子
    DF_oa = U_air*W*δ_air/A # U_air 混合区大气流速风速，cm/s,各变量需根据产地调查获得
    
    #表层土壤中污染物扩散进入室外空气的挥发因子，两个算法取较小值
    VF_suroa1 = rou_b/DF_oa*sqrt(4*Deff_ss*H/(pi*τ*3153600*K_sw*rou_b))*(1000)
    VF_suroa2 = d*rou_b/(DF_oa*τ*3156000)*(1000) # d为表层污染土壤层厚度，必须根据场地调查获得参数值;τ见G.1  
    VF_suroa = min(VF_suroa1, VF_suroa2)
    
    #敏感用地致癌效应
    IOVER_ca1 = VF_suroa*((DAIR_c*EFO_c*ED_c)/(BW_c*AT_ca)+(DAIR_a*EFO_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    IOVER_nc1 = VF_suroa*(DAIR_c*EFO_c*ED_c)/(BW_c*AT_nc)
    return(IOVER_ca1,IOVER_nc1)


# 暴露途径5 吸入室外空气中来自下层土壤的气态污染物
def IOVER2():
    #非饱和土层土壤中总孔隙体积比(F.2)、孔隙水体积比(F.3)、孔隙空气体积比(F.4)
    θ = 1-(ρ_b/ρ_s)   #ρ见G.1
    θ_ws = ρ_b*ρ_ws/ρ_w
    θ_as = θ-θ_ws 
    # 土壤中气态污染物的有效扩散系数(F.1)
    Deff_s = D_a*(θ_as**3.33)/(θ**2)+D_w*(θ_ws**3.33)/(H*θ**2) # D,H见B.2
    
    #土壤有机碳/土壤孔隙水分配系数
    K_d = K_oc*f_om/1700 #f_om为土壤有机质含量，需根据场地调查获得
    # 土壤-水污染物分配系数
    K_sw = (θ_ws+(K_d*ρ_b)+(H*θ_as))/ρ_b
    
    #室外空气中气态污染物扩散因子
    DF_oa = U_air*W*δ_air/A # U_air 混合区大气流速风速，cm/s,各变量需根据产地调查获得
    
    #下层土壤中污染物扩散进入室外空气的挥发因子
    VF_suboa1 = 1/((1+DF_oa*L_s/Deff_s)*K_sw/H)*1000
    VF_suboa2 = d_s*ρ_b/(DF_oa*τ*3153600)*1000
    VF_suboa = min(VF_suboa1, VF_suboa2)
    
    
    #敏感用地致癌效应
    IOVER_ca2 = VF_suboa*((DAIR_c*EFO_c*ED_c)/(BW_c*AT_ca)+(DAIR_a*EFO_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    IOVER_nc2 = VF_suboa*(DAIR_c*EFO_c*ED_c)/(BW_c*AT_nc)
    return(IOVER_ca2,IOVER_nc2)
    
    
# 暴露途径6 吸入室内空气中来自下层土壤的气态污染物
def IOVER3():
    #毛细管层中气态污染物的有效扩散系数
    Deff_cap = D_a*(D_acap**3.33)/(θ**2)+D_w*(D_wcap**3.33)/(H*θ**2)
    #气态污染物从地下水到表层土壤的有效扩散系数
    Deff_gws = L_gw/((h_cap/Deff_cap)+(h_v/Deff_s))
    #地下水中污染物扩散进入室外空气的挥发因子
    VF_gwoa = 1/((1+DF_oa*L_gw/Deff_gws)*1/H)*1000
    
    #室外空气中气态污染物扩散因子
    DF_oa = U_air*W*δ_air/A # U_air 混合区大气流速风速，cm/s,各变量需根据产地调查获得
    
    #敏感用地致癌效应
    IOVER_ca3 = VF_gwoa*((DAIR_c*EFO_c*ED_c)/(BW_c*AT_ca)+(DAIR_a*EFO_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    IOVER_nc3 = VF_gwoa*(DAIR_c*EFO_c*ED_c)/(BW_c*AT_nc)
    return(IOVER_ca3,IOVER_nc3)

    
# 暴露途径7 吸入室外空气中来自地下水的气态污染物
def IIVER1():
    #非饱和土层土壤中总孔隙体积比(F.2)、孔隙水体积比(F.3)、孔隙空气体积比(F.4)
    θ = 1-(ρ_b/ρ_s)   #ρ见G.1
    θ_ws = ρ_b*ρ_ws/ρ_w
    θ_as = θ-θ_ws 
    # 土壤中气态污染物的有效扩散系数(F.1)
    Deff_s = D_a*(θ_as**3.33)/(θ**2)+D_w*(θ_ws**3.33)/(H*θ**2) # D,H见B.2
    
    #流经地下室地板裂隙的对流空气流速
    R_crack = A_b*η/X_crack # η见G.1
    Q_s = 2*π*dP*K_v*X_crack/(μ_air*ln(2*Z_crack/R_crack)) #μ_air = 1.81*10^(-4),其余变量需测量
    if Q_s == 0:
        VF_subia1 = 1000/(K_sw/H*(1+Deff_s/(DF_ia*L_s)+Deff_s*L_crack/(Deff_s*L_s*η))*DF_ia/Deff_s*L_s)
    elif Q_s>0:
        ξ = Q_s*L_crack/(A_b*Deff_crack*η)
        VF_subia1 = 1000/(K_sw/H*(e**ξ+Deff_s/(DF_ia*L_s)+Deff_s*A_b/(Q_s*L_s)*(e**ξ-1))*DF_ia/Deff_s*L_s)
    VF_subia2 = d_s*ρ_b/(DF_ia*τ*31536000)*1000
    VF_subia = min(VF_subia1,VF_subia2)
    
    #敏感用地致癌效应
    IIVER_ca1 = VF_subia*((DAIR_c*EFI_c*ED_c)/(BW_c*AT_ca)+(DAIR_a*EFI_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    IIVER_nc1 = VF_subia*(DAIR_c*EFI_c*ED_c)/(BW_c*AT_nc)
    return(IIVER_ca1,IIVER_nc1)
    
    
# 暴露途径8 吸入室内空气中来自地下水的气态污染物
def IIVER2():
    #敏感用地致癌效应
    IIVER_ca2 = VF_gwia*((DAIR_c*EFI_c*ED_c)/(BW_c*AT_ca)+(DAIR_a*EFI_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    IIVER_nc2 = VF_gwia*(DAIR_c*EFI_c*ED_c)/(BW_c*AT_nc)
    return(IIVER_ca2,IIVER_nc2)


# 暴露途径9 饮用地下水
def CGWER():
    #敏感用地致癌效应
    CGWER_ca = ((GWCR_c*EF_c*ED_c)/(BW_c*AT_ca)+(GWCR_a*EF_a*ED_a)/(BW_a*AT_ca))
    #敏感用地非致癌效应
    CGWER_nc = (GWCR_c*EF_c*ED_c)/(BW_c*AT_nc)
    return(CGWER_ca,CGWER_nc)
    

a,b = OISER(OSIR_c,ED_c,EF_c,BW_c,OSIR_a,ED_a,EF_a,BW_a,ABS_0,AT_ca,AT_nc)