Image-based BRDF measurement
File_Hierarchy
-Calibration
	-some imgs for calibration
	-trackxxx.xml used to calibration camera track
	-in.xml used to calibration instrinsic parameters
	-pic_abced.xml input imgs list
	-default.xml to set up calibration software
-docs
	-ppt
	-report
-imgsRendered
	-used_to_find_BRDF.tga(high quality,1920*1080, used to generate final_data.txt)
	-pre_rendered(low quality,192*108, 200 pics, to show rotation)
-imgsTaken
	-Round1 used for calibration(with additional light)
	-Round2 used for BRDF measurement(with only point light)
-Software
REFERENCE:some code and standard from MIT_6.837 http://groups.csail.mit.edu/graphics/classes/6.837/F01/course-info.html
IN MAIN.CPP: Constants in 404/405/412 lines are used to do optimization
-final_data.txt(scatter point BRDF)

特别强调：
比展示时的优化：
有很多，真的有很多！首先我们做完了，不是fail的状态了！有了一组相对合理的数据。
数据为什么合理在报告的第六部分中有分析。
得到合理的数据所做的额外工作包括：
1.意识到两组照片（Round1中加入额外光源）更好calibrate。
2.于是调整棋盘，小车的位置重新拍了Round1和Round2的照片。
3.做了误差分析，用PhotoShop分析了pre-rendered的图和真实的图像素坐标误差差多少
4.并根据这个误差，在main.cpp软件中加了适当的优化方法（扩大像素搜索范围）对应的即为404、405、412行。

做出结果不易啊啊，希望老师、助教可以多给点分!!!超级感谢！！！

梁晨
2020.11.16
