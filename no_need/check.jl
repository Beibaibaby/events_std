using JLD2

# 加载 JLD2 文件
file_path = "/gpfs/data/doiron-lab/draco/weights_2024-03-20_12-31-30.jld2"
jldfile = jldopen(file_path, "r")  # 以只读模式打开

# 列出文件中所有的键
keyskk = keys(jldfile)

# 打印所有键
println("Keys stored in the JLD2 file:")
for key in keyskk
    println(key)
end

# 关闭 JLD2 文件
close(jldfile)
