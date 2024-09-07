function generateHeaderTemplate(filename, description)
% 生成 MATLAB 头注释模板
 
% 获取当前日期
dateStr = datestr(now, 'yyyy-mm-dd');
 
% 读取原始文件内容
content = fileread(filename);
 
% 创建头注释
header = sprintf('%%**************************************************************************\n');
header = [header, sprintf('%% 文件名: %s\n', filename)];
header = [header, sprintf('%% 版本: %s\n', "v1.0")];
header = [header, sprintf('%% 作者: %s\n', "Wuaoli-123")];
header = [header, sprintf('%% 联系方式: %s\n', "2713337051@qq.com")];
header = [header, sprintf('%% 日期: %s\n', dateStr)];
header = [header, sprintf('%% 描述: %s\n', description)];
header = [header, sprintf('%% 输入: %s\n', " ")];
header = [header, sprintf('%% 输出: %s\n', " ")];
header = [header, sprintf('%%**************************************************************************\n')];
 
% 将头注释和原始内容写入文件
fid = fopen(filename, 'w');
fprintf(fid, '%s\n', header);
fprintf(fid, '%s', content);
fclose(fid);
 
disp('头注释模板已生成。');
end