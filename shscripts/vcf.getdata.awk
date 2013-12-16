BEGIN { split(col,columns, ":"); split(formatname, formatfields, ":") }
$1 !~ /^#/ { \
   oidx = 1
   delete output
   split($9, format, ":")
   for (colidx=10; colidx<=NF; colidx++) {
       split($colidx, data, ":")
       for (i=1; i<=length(format); i++) {
       	   for (j=1; j<= length(formatfields); j++) {
       	       if (format[i] == formatfields[j]) {
	       	  output[oidx++] = data[i]
	       }  
	   }
       }
   }
   for (i=1; i<=length(columns); i++) {
       printf("%s\t",$(columns[i]))
   }
   for (i=1; i<=length(output); i++) {
       printf("%s\t",output[i])
   }
   printf("\n")
}

