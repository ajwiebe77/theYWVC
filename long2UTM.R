# O'Brien, J., 2012.
# https://stackoverflow.com/questions/9186496/determining-utm-zone-to-convert-from-longitude-latitude
long2UTM <- function(long) {
  (floor((long + 180)/6) %% 60) + 1
}
