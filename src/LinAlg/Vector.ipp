#pragma once

namespace LinAlg {
//******************************************************************************

template <typename T> Vector<T> Vector<T>::conj() const {
  static_assert(is_complex_v<T>, "conj() only available for complex Vector");
  std::vector<T> conj_data;
  conj_data.reserve(this->size());
  for (std::size_t i = 0; i < this->size(); ++i) {
    conj_data.push_back(std::conj(this->data()[i]));
  }
  return Vector<T>{std::move(conj_data)};
}

template <typename T> auto Vector<T>::real() const {
  static_assert(is_complex_v<T>, "real() only available for complex Vector");
  std::vector<typename T::value_type> real_data;
  real_data.reserve(this->size());
  for (std::size_t i = 0; i < this->size(); ++i) {
    real_data.push_back(std::real(this->data()[i]));
  }
  return Vector<typename T::value_type>{std::move(real_data)};
}
template <typename T> auto Vector<T>::imag() const {
  static_assert(is_complex_v<T>, "imag() only available for complex Vector");
  std::vector<typename T::value_type> imag_data;
  imag_data.reserve(this->size());
  for (std::size_t i = 0; i < this->size(); ++i) {
    imag_data.push_back(std::imag(this->data()[i]));
  }
  return Vector<typename T::value_type>{std::move(imag_data)};
}
template <typename T> auto Vector<T>::complex() const {
  static_assert(!is_complex_v<T>, "complex() only available for real Vector");
  // use move constructor to avoid default Matrix construction:
  std::vector<std::complex<T>> new_data;
  new_data.reserve(this->m_data.size());
  for (std::size_t i = 0; i < this->m_data.size(); ++i) {
    new_data.push_back(this->m_data[i]);
  }
  return Vector<std::complex<T>>{std::move(new_data)};
}

//******************************************************************************
template <typename T> Vector<T> &Vector<T>::operator+=(const Vector<T> &rhs) {
  assert(this->rows() == rhs.rows() && this->cols() == rhs.cols());
  using namespace qip::overloads;
  this->m_data += rhs.m_data;
  return *this;
}
template <typename T> Vector<T> &Vector<T>::operator-=(const Vector<T> rhs) {
  assert(this->rows() == rhs.rows() && this->cols() == rhs.cols());
  using namespace qip::overloads;
  this->m_data -= rhs.m_data;
  return *this;
}
template <typename T> Vector<T> &Vector<T>::operator*=(const T x) {
  using namespace qip::overloads;
  this->m_data *= x;
  return *this;
}
template <typename T> Vector<T> &Vector<T>::operator/=(const T x) {
  using namespace qip::overloads;
  this->m_data /= x;
  return *this;
}

//******************************************************************************
template <typename T> auto Vector<T>::as_gsl_view() {
  const auto size = std::max(this->rows(), this->cols());
  if constexpr (std::is_same_v<T, double>) {
    return gsl_vector_view_array(this->data(), size);
  } else if constexpr (std::is_same_v<T, float>) {
    return gsl_vector_float_view_array(this->data(), size);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    return gsl_vector_complex_view_array(
        reinterpret_cast<double *>(this->data()), size);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    return gsl_vector_complex_float_view_array(
        reinterpret_cast<float *>(this->data()), size);
  } else {
    assert(false &&
           "as_gsl_view only for double/float (or complex double/float)");
  }
}

template <typename T> auto Vector<T>::as_gsl_view() const {
  const auto size = std::max(this->rows(), this->cols());
  if constexpr (std::is_same_v<T, double>) {
    return gsl_vector_const_view_array(this->data(), size);
  } else if constexpr (std::is_same_v<T, float>) {
    return gsl_vector_float_const_view_array(this->data(), size);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    return gsl_vector_complex_const_view_array(
        reinterpret_cast<const double *>(this->data()), size);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    return gsl_vector_complex_float_const_view_array(
        reinterpret_cast<const float *>(this->data()), size);
  } else {
    assert(false &&
           "as_gsl_view only for double/float (or complex double/float)");
  }
}

} // namespace LinAlg
