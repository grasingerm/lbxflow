Get-ChildItem "..\ie_sims" -Filter *.js | `
Foreach-Object {
    $content = Get-Content $_.FullName

    #filter and save content to the original file
    $content | Where-Object {$_ -match 'step[49]'} | Set-Content $_.FullName
}
