set list = `\ls -d ?s_?`
echo $list
foreach file ($list)
   zip ${file}.zip $file/*
end
zip Alaska-PowellSources-v3.zip ?s_?.zip

