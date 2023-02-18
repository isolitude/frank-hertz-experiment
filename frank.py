from xml.dom.pulldom import END_ELEMENT
import random
import math
import numpy as np
import argon
import time


#常量定义
e=1.6*1e-19
m_e=9.1*1e-31
k=1.38*1e-23
T=293.15
p=1000
d_k2=0.02

#管内气体原子密度
n_gas=p/(k*T)

#模拟精度，最小时间单位
t_0=1e-10

#模拟电子数
n_e=1000

#定义类函数电子，包含速度，能量，碰撞截面等信息
class electron:

    v_e=0                        #粒子z轴速度 
    E_e=0                        #粒子能量
    X_e=0                        #粒子在k1k2极板间的位置
    section=0                    #总散射截面大小
    section1=0                   #激发碰撞截面11.83
    section2=0                   #弹性碰撞截面
    section3=0                   #激发碰撞截面14.09
    section4=0                   #激发碰撞截面14.30
    section5=0                   #电离碰撞15.76
    count=0.0                    #电离电子计数
    theta=0.0                    #电子速度方向，默认单位弧度


    def __init__(self,V_k1) -> None:
        #初始化速度分布，使用玻尔兹曼分布律，σ=（KT/m）^0.5,μ=0，的正太分布，对于复数的部分取绝对值
        sigma=math.sqrt(k*T/m_e)
        electron.v_e=math.sqrt(2*V_k1*e/m_e)+abs(np.random.normal(0,sigma,1))
        electron.X_e=0
        electron.count=0.0#电离碰撞计数
        electron.theta=0.0#初始速度方向

    def velocity(self,E):
        #更新速度
        delt_v=E*e*t_0/m_e #计算单位时间速度增量
        electron.X_e=electron.X_e+electron.v_e*t_0 #计算电子当前位置
        electron.v_e=electron.v_e+delt_v #计算电子当前速度
        #利用Vxy速度在电场下不变，列出守恒式，微分得(tanθ)dV+Vd(tanθ)=0
        #保证速度方向向z轴正向偏移
        electron.theta=electron.theta-delt_v*np.sin(electron.theta)*np.cos(electron.theta)/electron.v_e
    
    def cross(self):

        #更新能量
        electron.E_e=0.5*m_e*(electron.v_e/np.cos(electron.theta))**2

        #总散射截面
        electron.section=argon.totalsection(electron.E_e)

        #11.83eV
        #electron.section1=argon.excitationsection1(electron.E_e)

        #14.09eV
        #electron.section3=argon.excitationsection2(electron.E_e)

        #14.30eV
        #electron.section4=argon.excitationsection3(electron.E_e)

        #电离碰撞函数
        #electron.section5=argon.ionizationsection(electron.E_e)

    def impact(self,V_k2,V_a2):
        #调用此函数来判断是否发生碰撞
        #传入参数V_k2,V_a2用于判断电离电子，是否能到达电极被吸收
        #注意到除弹性碰撞外，都不会发生速度方向改变

        #碰撞概率的速度是总速度
        Poss=1-math.exp(-n_gas*electron.section*abs(electron.v_e/np.cos(electron.theta))*t_0)  #和管内气体有关

        #随机数判断碰撞
        if random.random()<=Poss:
            fate=random.random()
            electron.section1=argon.excitationsection1(electron.E_e)
            if fate<=(electron.section1/electron.section):
                
                #发生碰撞，电子能量被吸收11.83ev
                electron.v_e=math.sqrt(2*(electron.E_e-11.83*e)/m_e)*np.cos(electron.theta)  
                 
            else:
                electron.section3=argon.excitationsection2(electron.E_e)
                if fate<=((electron.section3+electron.section1)/electron.section):

                    #发生碰撞，电子能量吸收14.09eV
                    electron.v_e=math.sqrt(2*(electron.E_e-14.09*e)/m_e)*np.cos(electron.theta)  

                else:
                    electron.section4=argon.excitationsection3(electron.E_e)
                    if  fate<=((electron.section3+electron.section1+electron.section4)/electron.section):

                        #发生碰撞，电子能量吸收14.30eV
                        electron.v_e=math.sqrt(2*(electron.E_e-14.30*e)/m_e)*np.cos(electron.theta)  

                
                    else:
                        electron.section5=argon.ionizationsection(electron.E_e)
                        if fate<=((electron.section1+electron.section3+electron.section4+electron.section5)/electron.section):
                
                            #发生电离碰撞，导致出现一个电子，初速度为0，在当前位置，且不发生碰撞
                            electron.v_e=math.sqrt(2*(electron.E_e-15.76*e)/m_e)*np.cos(electron.theta)  


                            #判断当前电势，高于反向电压，则允许射出
                            #由于动量守恒，电离电子速度方向与原电子相反，平分剩余能量
                
                            if electron.X_e<=d_k2*(1-(V_a2-7.88+electron.E_e/e)/V_k2):
                                electron.count=electron.count+0.5
               
                

                        else:

                            #发生弹性碰撞
                            #加入随机角度分布
                            theta1=electron.theta
                            electron.theta=electron.theta+np.random.normal(0,0.1,1)
                            #更新z方向速度,用原速度除以原角度乘以现角度
                            electron.v_e=electron.v_e*np.cos(electron.theta)/np.cos(theta1)


#控制面板调用函数
def simulation(V_k2,V_a2,V_k1):
    n_e=1000
    v_2=math.sqrt(2*V_a2*e/m_e)
    E_k2=V_k2/d_k2
    n=0.0
    for i in range(0,n_e+1):
        el=electron(V_k1) #召唤一个电子
        while (el.X_e<d_k2)and(el.X_e>=0):
             #当电子运行到极板k2时跳出循环
             el.velocity(E_k2) #更新速度位置
             el.cross() #更新截面
             el.impact(V_k2,V_a2) #判断碰撞
        if el.v_e>=v_2:
            n=n+1
        n=n+el.count
        time.sleep(0.0001)
    return [V_k2,n]


def simulation1(args):
    V_k2,V_a2,V_k1=args[0],args[1],args[2]
    n_e=1000
    v_2=math.sqrt(2*V_a2*e/m_e)
    E_k2=V_k2/d_k2
    n=0.0
    for i in range(0,n_e+1):
        el=electron(V_k1) #召唤一个电子
        while (el.X_e<d_k2)and(el.X_e>=0):
             #当电子运行到极板k2时跳出循环
             el.velocity(E_k2) #更新速度位置
             el.cross() #更新截面
             el.impact(V_k2,V_a2) #判断碰撞
        if el.v_e>=v_2:
            n=n+1
        n=n+el.count
        
    print([V_k2,n])
    return [V_k2,n]


