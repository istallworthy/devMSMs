#' Title
#' This code comes from the aMNLFA package by Gottfredson and colleagues (2019)
#' @param df dataframe
#' @param outfile path to output file
#' @seealso [aMNLFA::aMNLFA.initial()]

#' @examples write.inp.file(df, outfile)
write.inp.file <- function(df, outfile) {

  #Remove final slash to allow for better use of file.path down the line.
  fixPath<-function(somecharacter) {
    lastcharacter<-base::substring(somecharacter,nchar(somecharacter))
    newcharacter<-base::substr(somecharacter,1,nchar(somecharacter)-1)
    outcharacter<-somecharacter
    if (lastcharacter=="/") {outcharacter<-newcharacter}
    outcharacter
  }

  #checking inputs
  if (!is.data.frame(df))
    stop("df argument needs to be a data.frame")
  if (!is.character(outfile))
    stop("outfile argument needs to be a character vector")
  # check if the outfile already exists, if it does clear it.
  if (file.exists(outfile)) file.remove(outfile)
  #replace missing values with ""
  df[is.na(df[,1]),1] <- ""
  #df[is.na(df)] = "" # RDS commented out above line and added this
  # next we iteratively write the .inp file from our input dataframe.
  for (i in seq(nrow(df))){
    # if number of characters is <= 90, write line and move to next loop
    if (nchar(df[i,1])<=90) {
      cat(df[i,1],
          sep = "\n",
          file = outfile,
          append = TRUE)
      # if characters between 90 and 160 we break it into two lines, at a space
      # note this can easily be extended to more than 160, if this solution works
    } else {
      # identify space closest to character limit 90 --this was originally only done for the first line --IS amended to iterate thru subsequent lines
      # first this locates all spaces in the current line

      #CREATING FIRST LINE OF 90 CHARACTERS OR FEWER
      sp <- as.matrix(stringr::str_locate_all(df[i,1], " ")[[1]])
      # then we select only spaces prior to the 90th character
      sp_resid=sp[sp[,1]>90,] #IS added to catch any excess to subsequently iterate on
      sp <- sp[sp[,1]<=90,1]
      # select the space closest to character limit
      sp <- sp[length(sp)]
      # cut vector at that space, and write as two lines
      cat(
        # write first line from 1 to the space "sp"
        substr(df[i,1], 1, sp),
        # write the second line from "sp" to the length of the vector
        # substr(df[i,1], sp, nchar(df[i,1])), #populating excess, assuming it is <90 chars
        sep = "\n",
        file = outfile,
        append = TRUE)

      #CREATING SECOND LINE OF 90 CHARACTERS OR FEWER (if necessary) added by IS
      #if the first chunking didn't leave less than 90 characters remaining for the second line, repeat
      if (max(sp_resid)-min(sp_resid)>90){ #if the remaining characters on this line are still longer than 90, do process again (added by IS)
        sp_resid2=sp_resid[sp_resid[,1]>(90+min(sp_resid)),] #IS added to catch any further excess
        sp2=sp_resid[sp_resid[,1]<(90+min(sp_resid)),1]
        sp2_end <- sp2[length(sp2)-1]
        sp2_start <- sp2[1]

        # cut vector at that space, and write as two lines
        cat(
          # write first line from start of new chunk to end of new chunk
          substr(df[i,1], sp, sp2_end),
          # substr(df[i,1], sp2_start, sp2_end),

          # write the second line from "sp" to the length of the vector
          # substr(df[i,1], sp2_end, nchar(df[i,1])),
          sep = "\n",
          file = outfile,
          append = TRUE)

        if (max(sp_resid2)-sp2_end>90){ #if the remaining characters on this line are still longer than 90, do process again (added by IS)
          sp_resid3=sp_resid2[sp_resid2[,1]>(90+min(sp_resid2)),] #IS added to catch any further excess
          sp3=sp_resid2[sp_resid2[,1]<(90+min(sp_resid2)),1]
          sp3_end <- sp3[length(sp3)-3]
          sp3_start <- sp3[1]

          # cut vector at that space, and write as two lines
          cat(
            # write first line from start of new chunk to end of new chunk
            substr(df[i,1], sp2_end, sp3_end),
            # write the second line from "sp" to the length of the vector
            substr(df[i,1], sp3_end, nchar(df[i,1])),
            sep = "\n",
            file = outfile,
            append = TRUE)
        } else {
          # cut vector at that space, and write as two lines
          cat(
            # write first line from start of new chunk to end of new chunk
            # substr(df[i,1], sp3_start, sp3_end),
            # write the second line from "sp" to the length of the vector
            substr(df[i,1], sp2_end, nchar(df[i,1])),
            sep = "\n",
            file = outfile,
            append = TRUE)
        }

      } else {
        cat(
          # write first line from 1 to the space "sp"
          # substr(df[i,1], 1, sp),
          # write the second line from "sp" to the length of the vector
          substr(df[i,1], sp, nchar(df[i,1])),
          sep = "\n",
          file = outfile,
          append = TRUE)
      }
    }
  }

}

fixPath<-function(somecharacter) {
  lastcharacter<-base::substring(somecharacter,nchar(somecharacter))
  newcharacter<-base::substr(somecharacter,1,nchar(somecharacter)-1)
  outcharacter<-somecharacter
  if (lastcharacter=="/") {outcharacter<-newcharacter}
  outcharacter
}
