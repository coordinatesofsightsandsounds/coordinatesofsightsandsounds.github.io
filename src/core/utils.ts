import { Coordinate } from "../components/app";

export const degreeToDecimalXY = ({
  degreeLat,
  minuteLat,
  secondLat,
  degreeLng,
  minuteLng,
  secondLng,
}: Coordinate): [number, number] => [
  degreeLng + minuteLng/60 + secondLng/3600,
  degreeLat + minuteLat/60 + secondLat/3600,
]

export const decimalToDegree = ([x, y]) => {
  const degreeLat = Math.floor(y)
  const minuteLat = Math.floor((y - degreeLat) * 60)
  const secondLat = Math.round((y - degreeLat - (minuteLat / 60)) * 36000) / 10
  const degreeLng = Math.floor(x)
  const minuteLng = Math.floor((x - degreeLng) * 60)
  const secondLng = Math.round((x - degreeLng - (minuteLng / 60)) * 36000) / 10

  return {
    degreeLat,
    minuteLat,
    secondLat,
    degreeLng,
    minuteLng,
    secondLng,
  }
}

globalThis.degreeToDecimalXY = degreeToDecimalXY
globalThis.decimalToDegree = decimalToDegree


export const degreeToString = ({
  degreeLat,
  minuteLat,
  secondLat,
  degreeLng,
  minuteLng,
  secondLng,
}: Coordinate) => `${Math.abs(degreeLat)}°${minuteLat}'${secondLat}"${degreeLat > 0 ? 'N' : 'S'} ${Math.abs(degreeLng)}°${minuteLng}'${secondLng}"${degreeLng > 0 ? 'E' : 'W'}`
