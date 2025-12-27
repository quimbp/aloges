!> @file module_nn.f90
!! @brief Simple fully-connected neural network (1 hidden layer) with optional normalization.
!!
!! @details
!! This module implements a small feed-forward neural network intended for
!! function approximation / regression on small datasets.
!!
!! Features:
!! - 1 hidden layer fully-connected network.
!! - Selectable hidden and output activation functions.
!! - Optional L2 regularization on weights.
!! - Optional z-score normalization (standardization) for inputs/targets using
!!   statistics computed from the training set.
!! - Training via L-BFGS-B (external routine `setulb`).
!!
!! Data layout conventions:
!! - Inputs are stored as a matrix of shape (n_inputs, n_samples).
!! - Targets are stored as a matrix of shape (n_outputs, n_samples).
!! - Samples are column-major: sample i is `inputs(:,i)`.
!!
!! Normalization conventions:
!! - `*_original` arrays keep original (unnormalized) data for reference.
!! - `training_inputs`, `training_targets`, etc. are *working copies* that may be
!!   normalized and used for training.
!! - `predict()` ALWAYS accepts inputs in ORIGINAL space and ALWAYS returns
!!   outputs in ORIGINAL space.
!!
!! @warning
!! Any routine that evaluates predictions vs. ORIGINAL targets must use
!! `predict()` (not `forward()` directly), unless the caller explicitly provides
!! already-normalized inputs/targets and knows what they are doing.
!!
!! @note
!! This module depends on `aloges module_type` for the `dp` kind and the 
!! external L-BFGS-B `setulb` routine (or interface).

module module_nn

use module_types, only: dp
implicit none

!> Activation function ID: hyperbolic tangent.
integer, parameter :: ACTIVATION_TANH    = 1
!> Activation function ID: logistic sigmoid.
integer, parameter :: ACTIVATION_SIGMOID = 2
!> Activation function ID: ReLU.
integer, parameter :: ACTIVATION_RELU    = 3
!> Activation function ID: linear (identity).
integer, parameter :: ACTIVATION_LINEAR  = 4

!> @brief Neural network object (single hidden layer).
!!
!! @details
!! Stores architecture, parameters, datasets, and normalization statistics.
!!
!! Intended usage pattern:
!! 1. `call net%init(n_in, n_hidden, n_out[, seed])`
!! 2. `call net%set_training_data(inputs, targets)`
!! 3. Optional: `call net%set_validation_data(val_inputs, val_targets)`
!! 4. Optional: `call net%normalize_data()`  (recommended for many problems)
!! 5. `call net%train([tol], [max_iter], [verbose])`
!! 6. Inference: `y = net%predict(x)`
!!
!! @note
!! When normalization is enabled (`is_normalized=.true.`), training optimizes a
!! loss in NORMALIZED space (because `training_inputs/targets` are normalized).
!! `predict()` transparently maps between original and normalized spaces.
type type_nn
  !> @name Architecture
  !! @{
  integer :: n_inputs   !< Number of input features.
  integer :: n_hidden   !< Number of hidden units.
  integer :: n_outputs  !< Number of outputs.
  integer :: n_params   !< Total number of trainable parameters.
  !! @}

  !> @name Parameters (weights/biases)
  !! @{
  real(dp), allocatable :: weights_ih(:,:)  !< Input-to-hidden weights, shape (n_inputs, n_hidden).
  real(dp), allocatable :: weights_ho(:,:)  !< Hidden-to-output weights, shape (n_hidden, n_outputs).
  real(dp), allocatable :: bias_h(:)        !< Hidden biases, length n_hidden.
  real(dp), allocatable :: bias_o(:)        !< Output biases, length n_outputs.
  !! @}

  !> @name Working buffers
  !! @{
  real(dp), allocatable :: hidden(:)  !< Hidden activations for last forward pass.
  real(dp), allocatable :: output(:)  !< Output activations for last forward pass.
  !! @}

  !> @name Activation configuration
  !! @{
  integer :: hidden_activation = ACTIVATION_TANH   !< Activation for hidden layer.
  integer :: output_activation = ACTIVATION_LINEAR !< Activation for output layer.
  !! @}

  !> @name Regularization
  !! @{
  real(dp) :: l2_lambda = 0.0_dp  !< L2 regularization strength (>=0).
  !! @}

  !> @name Training dataset (working copies; may be normalized)
  !! @{
  integer :: n_samples  !< Number of training samples.
  real(dp), allocatable :: training_inputs(:,:)   !< Training inputs, shape (n_inputs, n_samples).
  real(dp), allocatable :: training_targets(:,:)  !< Training targets, shape (n_outputs, n_samples).
  !! @}

  !> @name Training dataset (original, unnormalized)
  !! @details
  !! These are preserved for reference and to compute normalization statistics.
  !! @{
  real(dp), allocatable :: training_inputs_original(:,:)   !< Original training inputs.
  real(dp), allocatable :: training_targets_original(:,:)  !< Original training targets.
  !! @}

  !> @name Validation dataset (working copies; may be normalized)
  !! @{
  integer :: n_validation  !< Number of validation samples.
  real(dp), allocatable :: validation_inputs(:,:)   !< Validation inputs, shape (n_inputs, n_validation).
  real(dp), allocatable :: validation_targets(:,:)  !< Validation targets, shape (n_outputs, n_validation).
  !! @}

  !> @name Validation dataset (original, unnormalized)
  !! @{
  real(dp), allocatable :: validation_inputs_original(:,:)   !< Original validation inputs.
  real(dp), allocatable :: validation_targets_original(:,:)  !< Original validation targets.
  !! @}

  !> @name Training history (currently unused)
  !! @details
  !! Present for future extensions (epoch-based training).
  !! @{
  integer :: n_epochs_trained = 0
  real(dp), allocatable :: train_error_history(:)
  real(dp), allocatable :: val_error_history(:)
  !! @}

  !> @name Normalization parameters
  !! @details
  !! If `is_normalized=.true.`, inputs are standardized via:
  !! $$x_{norm} = (x - mean)/std$$
  !! and outputs are de-standardized via:
  !! $$y = y_{norm} * std + mean$$
  !! @{
  logical :: is_normalized = .false.         !< Whether normalization stats are available and enabled.
  real(dp), allocatable :: input_mean(:)     !< Mean for each input dimension.
  real(dp), allocatable :: input_std(:)      !< Std for each input dimension (>=epsilon).
  real(dp), allocatable :: output_mean(:)    !< Mean for each output dimension.
  real(dp), allocatable :: output_std(:)     !< Std for each output dimension (>=epsilon).
  !! @}

contains
  !> @name Public type-bound procedures
  !! @{
  procedure :: init              => nn_initialize
  procedure :: forward           => nn_forward
  procedure :: destroy           => nn_destroy
  procedure :: set_training_data => nn_set_training_data
  procedure :: set_validation_data => nn_set_validation_data
  procedure :: vector_to_network => nn_vector_to_network
  procedure :: network_to_vector => nn_network_to_vector
  procedure :: predict           => nn_predict
  procedure :: evaluate          => nn_evaluate
  procedure :: compute_metrics   => nn_compute_metrics
  procedure :: normalize_data    => nn_normalize_data
  procedure :: denormalize_output => nn_denormalize_output
  procedure :: save_model        => nn_save_model
  procedure :: load_model        => nn_load_model
  procedure :: set_activation    => nn_set_activation
  procedure :: set_regularization => nn_set_regularization
  procedure :: train             => nn_train_lbfgsb
  procedure :: check_gradient    => nn_check_gradient
  !! @}
end type type_nn

!> @brief Pointer used to guard against simultaneous training callbacks.
!! @details
!! The L-BFGS-B routine (`setulb`) often uses callbacks. This module uses a
!! single global pointer to the currently-optimized network, so concurrent
!! training is prevented.
type(type_nn), pointer :: current_nn => null()

!> @brief Guard flag to prevent multiple networks being trained simultaneously.
logical :: callback_active = .false.

contains

  !> Initialize the neural network and allocate parameters/buffers.
  !!
  !! @param[inout] this  Network object to initialize.
  !! @param[in]    n_in  Number of inputs (features).
  !! @param[in]    n_hid Number of hidden units.
  !! @param[in]    n_out Number of outputs.
  !! @param[in]    seed  Optional RNG seed for reproducible initialization.
  !!
  !! @details
  !! Allocates weights/biases and initializes weights with a Xavier/Glorot-like
  !! uniform distribution. Biases are initialized to zero.
  !!
  !! @warning
  !! This routine does not clear datasets automatically; if you reuse the same
  !! instance, call `destroy()` if you want to free all memory.
  subroutine nn_initialize(this, n_in, n_hid, n_out, seed)
    class(type_nn), intent(inout) :: this
    integer, intent(in) :: n_in, n_hid, n_out
    integer, intent(in), optional :: seed
    integer :: seed_size, i
    integer, allocatable :: seed_array(:)

    if (n_in <= 0 .or. n_hid <= 0 .or. n_out <= 0) then
      error stop "Error: All network dimensions must be positive"
    end if

    this%n_inputs  = n_in
    this%n_hidden  = n_hid
    this%n_outputs = n_out
    this%n_params  = n_in * n_hid + n_hid * n_out + n_hid + n_out

    if (present(seed)) then
      call random_seed(size=seed_size)
      allocate(seed_array(seed_size))
      do i = 1, seed_size
        seed_array(i) = seed + i - 1
      end do
      call random_seed(put=seed_array)
      deallocate(seed_array)
    end if

    if (allocated(this%weights_ih)) deallocate(this%weights_ih)
    if (allocated(this%weights_ho)) deallocate(this%weights_ho)
    if (allocated(this%bias_h))     deallocate(this%bias_h)
    if (allocated(this%bias_o))     deallocate(this%bias_o)
    if (allocated(this%hidden))     deallocate(this%hidden)
    if (allocated(this%output))     deallocate(this%output)

    allocate(this%weights_ih(n_in, n_hid))
    allocate(this%weights_ho(n_hid, n_out))
    allocate(this%bias_h(n_hid))
    allocate(this%bias_o(n_out))
    allocate(this%hidden(n_hid))
    allocate(this%output(n_out))

    call random_number(this%weights_ih)
    call random_number(this%weights_ho)
    call random_number(this%bias_h)
    call random_number(this%bias_o)

    this%weights_ih = (this%weights_ih - 0.5_dp) * 2.0_dp * sqrt(6.0_dp / real(n_in + n_hid, dp))
    this%weights_ho = (this%weights_ho - 0.5_dp) * 2.0_dp * sqrt(6.0_dp / real(n_hid + n_out, dp))
    this%bias_h = 0.0_dp
    this%bias_o = 0.0_dp
  end subroutine nn_initialize

  !> Set activation functions for hidden and output layers.
  !!
  !! @param[inout] this        Network object.
  !! @param[in]    hidden_act  Activation ID for hidden layer.
  !! @param[in]    output_act  Activation ID for output layer.
  !!
  !! @details
  !! Valid activation IDs: `ACTIVATION_TANH`, `ACTIVATION_SIGMOID`,
  !! `ACTIVATION_RELU`, `ACTIVATION_LINEAR`.
  subroutine nn_set_activation(this, hidden_act, output_act)
    class(type_nn), intent(inout) :: this
    integer, intent(in) :: hidden_act, output_act

    if (hidden_act < 1 .or. hidden_act > 4) then
      error stop "Invalid hidden activation function"
    end if
    if (output_act < 1 .or. output_act > 4) then
      error stop "Invalid output activation function"
    end if

    this%hidden_activation = hidden_act
    this%output_activation = output_act
  end subroutine nn_set_activation

  !> Set L2 regularization strength.
  !!
  !! @param[inout] this      Network object.
  !! @param[in]    l2_lambda Non-negative regularization coefficient.
  !!
  !! @details
  !! Regularization is applied only to weights (not biases) in
  !! `compute_loss_and_gradient`.
  subroutine nn_set_regularization(this, l2_lambda)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: l2_lambda

    if (l2_lambda < 0.0_dp) then
      error stop "L2 regularization parameter must be non-negative"
    end if
    this%l2_lambda = l2_lambda
  end subroutine nn_set_regularization

  !> Apply activation function (scalar).
  !!
  !! @param[in] x               Pre-activation scalar.
  !! @param[in] activation_type Activation ID.
  !! @return y                  Activated scalar.
  !!
  !! @note
  !! This is a scalar helper used inside `nn_forward`.
  function apply_activation(x, activation_type) result(y)
    real(dp), intent(in) :: x
    integer, intent(in) :: activation_type
    real(dp) :: y

    select case(activation_type)
      case(ACTIVATION_TANH)
        y = tanh(x)
      case(ACTIVATION_SIGMOID)
        y = 1.0_dp / (1.0_dp + exp(-x))
      case(ACTIVATION_RELU)
        y = max(0.0_dp, x)
      case(ACTIVATION_LINEAR)
        y = x
      case default
        y = x
    end select
  end function apply_activation

  !> Set training data (stores both original and working copies).
  !!
  !! @param[inout] this    Network object.
  !! @param[in]    inputs  Training inputs (n_inputs, n_samples) in ORIGINAL space.
  !! @param[in]    targets Training targets (n_outputs, n_samples) in ORIGINAL space.
  !!
  !! @details
  !! This routine:
  !! - allocates `*_original` arrays and copies the provided data,
  !! - allocates working copies (`training_inputs`, `training_targets`) and
  !!   initially sets them equal to original data.
  !!
  !! If `normalize_data()` is later called, the working copies will be replaced
  !! with normalized values, but the original arrays remain unchanged.
  subroutine nn_set_training_data(this, inputs, targets)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:,:), targets(:,:)

    if (size(inputs, 1) /= this%n_inputs) then
      error stop "Error: Input dimensions do not match network architecture"
    end if
    if (size(targets, 1) /= this%n_outputs) then
      error stop "Error: Target dimensions do not match network architecture"
    end if
    if (size(inputs, 2) /= size(targets, 2)) then
      error stop "Error: Number of samples must match in inputs and targets"
    end if
    if (size(inputs, 2) == 0) then
      error stop "Error: No training samples provided"
    end if

    this%n_samples = size(inputs, 2)

    if (allocated(this%training_inputs))          deallocate(this%training_inputs)
    if (allocated(this%training_targets))         deallocate(this%training_targets)
    if (allocated(this%training_inputs_original)) deallocate(this%training_inputs_original)
    if (allocated(this%training_targets_original))deallocate(this%training_targets_original)

    allocate(this%training_inputs(this%n_inputs, this%n_samples))
    allocate(this%training_targets(this%n_outputs, this%n_samples))
    allocate(this%training_inputs_original(this%n_inputs, this%n_samples))
    allocate(this%training_targets_original(this%n_outputs, this%n_samples))

    this%training_inputs_original = inputs
    this%training_targets_original = targets

    this%training_inputs = inputs
    this%training_targets = targets
  end subroutine nn_set_training_data

  !> Set validation data (stores both original and working copies).
  !!
  !! @param[inout] this    Network object.
  !! @param[in]    inputs  Validation inputs (n_inputs, n_validation) in ORIGINAL space.
  !! @param[in]    targets Validation targets (n_outputs, n_validation) in ORIGINAL space.
  !!
  !! @details
  !! If `normalize_data()` is called, validation working copies are normalized
  !! using training-set statistics.
  subroutine nn_set_validation_data(this, inputs, targets)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:,:), targets(:,:)

    if (size(inputs, 1) /= this%n_inputs) then
      error stop "Error: Validation input dimensions do not match network"
    end if
    if (size(targets, 1) /= this%n_outputs) then
      error stop "Error: Validation target dimensions do not match network"
    end if
    if (size(inputs, 2) /= size(targets, 2)) then
      error stop "Error: Number of validation samples must match"
    end if

    this%n_validation = size(inputs, 2)

    if (allocated(this%validation_inputs))          deallocate(this%validation_inputs)
    if (allocated(this%validation_targets))         deallocate(this%validation_targets)
    if (allocated(this%validation_inputs_original)) deallocate(this%validation_inputs_original)
    if (allocated(this%validation_targets_original))deallocate(this%validation_targets_original)

    allocate(this%validation_inputs(this%n_inputs, this%n_validation))
    allocate(this%validation_targets(this%n_outputs, this%n_validation))
    allocate(this%validation_inputs_original(this%n_inputs, this%n_validation))
    allocate(this%validation_targets_original(this%n_outputs, this%n_validation))

    this%validation_inputs_original = inputs
    this%validation_targets_original = targets

    this%validation_inputs = inputs
    this%validation_targets = targets
  end subroutine nn_set_validation_data

  !> Compute normalization statistics on training set and normalize working copies.
  !!
  !! @param[inout] this Network object.
  !!
  !! @details
  !! Computes per-feature mean/std from `training_inputs_original` and per-output
  !! mean/std from `training_targets_original`.
  !!
  !! Produces normalized working copies:
  !! - `training_inputs`, `training_targets` become normalized
  !! - validation working copies (if present) are normalized using TRAINING stats
  !!
  !! @note
  !! Standard deviation is clamped to at least `epsilon` to avoid division by zero.
  !!
  !! @warning
  !! You should call `normalize_data()` after setting training (and optionally
  !! validation) data and before training, if you want normalized training.
  subroutine nn_normalize_data(this)
    class(type_nn), intent(inout) :: this
    integer :: i
    real(dp) :: epsilon = 1.0e-8_dp

    if (.not. allocated(this%training_inputs_original)) then
      error stop "Error: Training data must be set before normalization"
    end if

    if (allocated(this%input_mean))  deallocate(this%input_mean)
    if (allocated(this%input_std))   deallocate(this%input_std)
    if (allocated(this%output_mean)) deallocate(this%output_mean)
    if (allocated(this%output_std))  deallocate(this%output_std)

    allocate(this%input_mean(this%n_inputs))
    allocate(this%input_std(this%n_inputs))
    allocate(this%output_mean(this%n_outputs))
    allocate(this%output_std(this%n_outputs))

    do i = 1, this%n_inputs
      this%input_mean(i) = sum(this%training_inputs_original(i,:)) / this%n_samples
      this%input_std(i)  = sqrt(sum((this%training_inputs_original(i,:) - this%input_mean(i))**2) / this%n_samples)
      if (this%input_std(i) < epsilon) this%input_std(i) = 1.0_dp
      this%training_inputs(i,:) = (this%training_inputs_original(i,:) - this%input_mean(i)) / this%input_std(i)
    end do

    do i = 1, this%n_outputs
      this%output_mean(i) = sum(this%training_targets_original(i,:)) / this%n_samples
      this%output_std(i)  = sqrt(sum((this%training_targets_original(i,:) - this%output_mean(i))**2) / this%n_samples)
      if (this%output_std(i) < epsilon) this%output_std(i) = 1.0_dp
      this%training_targets(i,:) = (this%training_targets_original(i,:) - this%output_mean(i)) / this%output_std(i)
    end do

    if (allocated(this%validation_inputs_original)) then
      do i = 1, this%n_inputs
        this%validation_inputs(i,:) = (this%validation_inputs_original(i,:) - this%input_mean(i)) / this%input_std(i)
      end do
      do i = 1, this%n_outputs
        this%validation_targets(i,:) = (this%validation_targets_original(i,:) - this%output_mean(i)) / this%output_std(i)
      end do
    end if

    this%is_normalized = .true.
  end subroutine nn_normalize_data

  !> Convert normalized outputs back to original output space.
  !!
  !! @param[in] this              Network object.
  !! @param[in] normalized_output Output vector in normalized space.
  !! @return output               Output vector in original space.
  !!
  !! @details
  !! If `is_normalized=.false.`, this is the identity mapping.
  function nn_denormalize_output(this, normalized_output) result(output)
    class(type_nn), intent(in) :: this
    real(dp), intent(in) :: normalized_output(:)
    real(dp) :: output(size(normalized_output))
    integer :: i

    if (.not. this%is_normalized) then
      output = normalized_output
      return
    end if

    do i = 1, size(normalized_output)
      output(i) = normalized_output(i) * this%output_std(i) + this%output_mean(i)
    end do
  end function nn_denormalize_output

  !> Perform a forward pass using the provided input vector.
  !!
  !! @param[inout] this   Network object.
  !! @param[in]    inputs Input vector of length `n_inputs`.
  !!
  !! @details
  !! Writes results into internal buffers:
  !! - `this%hidden`
  !! - `this%output`
  !!
  !! @warning
  !! This routine does NOT apply normalization. If the network has been
  !! normalized, the caller must supply normalized inputs. Prefer `predict()`
  !! for external inference in original space.
  subroutine nn_forward(this, inputs)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:)
    integer :: i
    real(dp) :: z

    if (size(inputs) /= this%n_inputs) then
      error stop "Error: Input size does not match network architecture"
    end if

    do i = 1, this%n_hidden
      z = dot_product(inputs, this%weights_ih(:,i)) + this%bias_h(i)
      this%hidden(i) = apply_activation(z, this%hidden_activation)
    end do

    do i = 1, this%n_outputs
      z = dot_product(this%hidden, this%weights_ho(:,i)) + this%bias_o(i)
      this%output(i) = apply_activation(z, this%output_activation)
    end do
  end subroutine nn_forward

  !> Serialize network parameters into a flat vector.
  !!
  !! @param[in]  this   Network object.
  !! @param[out] params Flat parameter vector of length `n_params`.
  !!
  !! @details
  !! Layout is:
  !! 1. `weights_ih` flattened
  !! 2. `weights_ho` flattened
  !! 3. `bias_h`
  !! 4. `bias_o`
  subroutine nn_network_to_vector(this, params)
    class(type_nn), intent(in) :: this
    real(dp), intent(out) :: params(:)
    integer :: idx

    if (size(params) /= this%n_params) then
      error stop "Error: Parameter vector size mismatch"
    end if

    idx = 1
    params(idx:idx+this%n_inputs*this%n_hidden-1) = reshape(this%weights_ih, [this%n_inputs*this%n_hidden])
    idx = idx + this%n_inputs*this%n_hidden

    params(idx:idx+this%n_hidden*this%n_outputs-1) = reshape(this%weights_ho, [this%n_hidden*this%n_outputs])
    idx = idx + this%n_hidden*this%n_outputs

    params(idx:idx+this%n_hidden-1) = this%bias_h
    idx = idx + this%n_hidden

    params(idx:idx+this%n_outputs-1) = this%bias_o
  end subroutine nn_network_to_vector

  !> Load network parameters from a flat vector.
  !!
  !! @param[inout] this   Network object.
  !! @param[in]    params Flat parameter vector of length `n_params`.
  subroutine nn_vector_to_network(this, params)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: params(:)
    integer :: idx

    if (size(params) /= this%n_params) then
      error stop "Error: Parameter vector size mismatch"
    end if

    idx = 1
    this%weights_ih = reshape(params(idx:idx+this%n_inputs*this%n_hidden-1), [this%n_inputs, this%n_hidden])
    idx = idx + this%n_inputs*this%n_hidden

    this%weights_ho = reshape(params(idx:idx+this%n_hidden*this%n_outputs-1), [this%n_hidden, this%n_outputs])
    idx = idx + this%n_hidden*this%n_outputs

    this%bias_h = params(idx:idx+this%n_hidden-1)
    idx = idx + this%n_hidden

    this%bias_o = params(idx:idx+this%n_outputs-1)
  end subroutine nn_vector_to_network

  !> Predict outputs for a single input vector in ORIGINAL space.
  !!
  !! @param[inout] this   Network object.
  !! @param[in]    inputs Input vector (length n_inputs) in ORIGINAL space.
  !! @return prediction   Output vector (length n_outputs) in ORIGINAL space.
  !!
  !! @details
  !! If `is_normalized=.true.`:
  !! - inputs are normalized using training-set statistics
  !! - forward pass is executed in normalized space
  !! - outputs are de-normalized before returning
  function nn_predict(this, inputs) result(prediction)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:)
    real(dp) :: prediction(this%n_outputs)
    real(dp), allocatable :: normalized_inputs(:)
    integer :: i

    if (size(inputs) /= this%n_inputs) then
      error stop "Error: predict(): input dimension mismatch"
    end if

    if (this%is_normalized) then
      allocate(normalized_inputs(this%n_inputs))
      do i = 1, this%n_inputs
        normalized_inputs(i) = (inputs(i) - this%input_mean(i)) / this%input_std(i)
      end do
      call this%forward(normalized_inputs)
      prediction = this%denormalize_output(this%output)
      deallocate(normalized_inputs)
    else
      call this%forward(inputs)
      prediction = this%output
    end if
  end function nn_predict

  !> Evaluate RMSE on a dataset in ORIGINAL space.
  !!
  !! @param[inout] this    Network object.
  !! @param[in]    inputs  Input matrix (n_inputs, n_samples) in ORIGINAL space.
  !! @param[in]    targets Target matrix (n_outputs, n_samples) in ORIGINAL space.
  !! @return rmse          Scalar RMSE.
  !!
  !! @details
  !! Uses `predict()` to ensure correct handling of normalization.
  function nn_evaluate(this, inputs, targets) result(rmse)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:,:), targets(:,:)
    real(dp) :: rmse
    integer :: i, n_samples
    real(dp) :: error_sum
    real(dp) :: pred(this%n_outputs)

    if (size(inputs,1) /= this%n_inputs) then
      error stop "Error: evaluate(): input dimension mismatch"
    end if
    if (size(targets,1) /= this%n_outputs) then
      error stop "Error: evaluate(): target dimension mismatch"
    end if
    if (size(inputs,2) /= size(targets,2)) then
      error stop "Error: evaluate(): sample count mismatch"
    end if

    n_samples = size(inputs, 2)
    if (n_samples <= 0) then
      error stop "Error: evaluate(): no samples"
    end if

    error_sum = 0.0_dp
    do i = 1, n_samples
      pred = this%predict(inputs(:,i))
      error_sum = error_sum + sum((targets(:,i) - pred)**2)
    end do

    rmse = sqrt(error_sum / (n_samples * this%n_outputs))
  end function nn_evaluate

  !> Compute regression metrics (RMSE, MAE, R^2) in ORIGINAL space.
  !!
  !! @param[inout] this    Network object.
  !! @param[in]    inputs  Input matrix (n_inputs, n_samples) in ORIGINAL space.
  !! @param[in]    targets Target matrix (n_outputs, n_samples) in ORIGINAL space.
  !! @param[out]   rmse    Root-mean-square error.
  !! @param[out]   mae     Mean absolute error.
  !! @param[out]   r2      Coefficient of determination (R^2).
  !!
  !! @details
  !! Uses `predict()` so that normalization is handled correctly.
  subroutine nn_compute_metrics(this, inputs, targets, rmse, mae, r2)
    class(type_nn), intent(inout) :: this
    real(dp), intent(in) :: inputs(:,:), targets(:,:)
    real(dp), intent(out) :: rmse, mae, r2
    integer :: i, n_samples
    real(dp) :: error_sum, abs_error_sum, ss_tot, ss_res
    real(dp) :: target_mean(this%n_outputs)
    real(dp) :: pred(this%n_outputs)

    if (size(inputs,1) /= this%n_inputs) then
      error stop "Error: compute_metrics(): input dimension mismatch"
    end if
    if (size(targets,1) /= this%n_outputs) then
      error stop "Error: compute_metrics(): target dimension mismatch"
    end if
    if (size(inputs,2) /= size(targets,2)) then
      error stop "Error: compute_metrics(): sample count mismatch"
    end if

    n_samples = size(inputs, 2)
    error_sum = 0.0_dp
    abs_error_sum = 0.0_dp
    ss_res = 0.0_dp
    ss_tot = 0.0_dp

    do i = 1, this%n_outputs
      target_mean(i) = sum(targets(i,:)) / n_samples
    end do

    do i = 1, n_samples
      pred = this%predict(inputs(:,i))
      error_sum = error_sum + sum((targets(:,i) - pred)**2)
      abs_error_sum = abs_error_sum + sum(abs(targets(:,i) - pred))
      ss_res = ss_res + sum((targets(:,i) - pred)**2)
      ss_tot = ss_tot + sum((targets(:,i) - target_mean)**2)
    end do

    rmse = sqrt(error_sum / (n_samples * this%n_outputs))
    mae  = abs_error_sum / (n_samples * this%n_outputs)

    if (ss_tot > 0.0_dp) then
      r2 = 1.0_dp - ss_res / ss_tot
    else
      r2 = 0.0_dp
    end if
  end subroutine nn_compute_metrics

  !> Save network architecture, parameters, and (optional) normalization stats.
  !!
  !! @param[in] this     Network object.
  !! @param[in] filename Output file path.
  !!
  !! @details
  !! Uses unformatted Fortran binary output. The file is not portable across
  !! compilers/platforms unless you control those details.
  subroutine nn_save_model(this, filename)
    class(type_nn), intent(in) :: this
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    open(newunit=unit, file=trim(filename), form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      error stop "Error: Cannot open file for writing"
    end if

    write(unit) this%n_inputs, this%n_hidden, this%n_outputs
    write(unit) this%hidden_activation, this%output_activation
    write(unit) this%l2_lambda
    write(unit) this%weights_ih
    write(unit) this%weights_ho
    write(unit) this%bias_h
    write(unit) this%bias_o

    if (this%is_normalized) then
      write(unit) .true.
      write(unit) this%input_mean, this%input_std
      write(unit) this%output_mean, this%output_std
    else
      write(unit) .false.
    end if

    close(unit)
  end subroutine nn_save_model

  !> Load network architecture, parameters, and (optional) normalization stats.
  !!
  !! @param[inout] this     Network object to populate.
  !! @param[in]    filename File to read.
  !!
  !! @warning
  !! This routine calls `init()` internally to allocate arrays and reset the
  !! architecture based on the file contents.
  subroutine nn_load_model(this, filename)
    class(type_nn), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer :: unit, ios, n_in, n_hid, n_out
    logical :: normalized

    open(newunit=unit, file=trim(filename), form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      error stop "Error: Cannot open file for reading"
    end if

    read(unit) n_in, n_hid, n_out
    call this%init(n_in, n_hid, n_out)

    read(unit) this%hidden_activation, this%output_activation
    read(unit) this%l2_lambda
    read(unit) this%weights_ih
    read(unit) this%weights_ho
    read(unit) this%bias_h
    read(unit) this%bias_o

    read(unit) normalized
    if (normalized) then
      if (allocated(this%input_mean))  deallocate(this%input_mean)
      if (allocated(this%input_std))   deallocate(this%input_std)
      if (allocated(this%output_mean)) deallocate(this%output_mean)
      if (allocated(this%output_std))  deallocate(this%output_std)
      allocate(this%input_mean(n_in), this%input_std(n_in))
      allocate(this%output_mean(n_out), this%output_std(n_out))
      read(unit) this%input_mean, this%input_std
      read(unit) this%output_mean, this%output_std
      this%is_normalized = .true.
    else
      this%is_normalized = .false.
    end if

    close(unit)
  end subroutine nn_load_model

  !> Deallocate all allocatable components of the network.
  !!
  !! @param[inout] this Network object.
  !!
  !! @details
  !! Safe to call multiple times.
  subroutine nn_destroy(this)
    class(type_nn), intent(inout) :: this

    if (allocated(this%weights_ih))             deallocate(this%weights_ih)
    if (allocated(this%weights_ho))             deallocate(this%weights_ho)
    if (allocated(this%bias_h))                 deallocate(this%bias_h)
    if (allocated(this%bias_o))                 deallocate(this%bias_o)
    if (allocated(this%hidden))                 deallocate(this%hidden)
    if (allocated(this%output))                 deallocate(this%output)
    if (allocated(this%training_inputs))        deallocate(this%training_inputs)
    if (allocated(this%training_targets))       deallocate(this%training_targets)
    if (allocated(this%training_inputs_original)) deallocate(this%training_inputs_original)
    if (allocated(this%training_targets_original))deallocate(this%training_targets_original)
    if (allocated(this%validation_inputs))      deallocate(this%validation_inputs)
    if (allocated(this%validation_targets))     deallocate(this%validation_targets)
    if (allocated(this%validation_inputs_original)) deallocate(this%validation_inputs_original)
    if (allocated(this%validation_targets_original))deallocate(this%validation_targets_original)
    if (allocated(this%train_error_history))    deallocate(this%train_error_history)
    if (allocated(this%val_error_history))      deallocate(this%val_error_history)
    if (allocated(this%input_mean))             deallocate(this%input_mean)
    if (allocated(this%input_std))              deallocate(this%input_std)
    if (allocated(this%output_mean))            deallocate(this%output_mean)
    if (allocated(this%output_std))             deallocate(this%output_std)
  end subroutine nn_destroy

  !> Train using the L-BFGS-B optimizer (external `setulb`).
  !!
  !! @param[inout] this     Network object.
  !! @param[in]    tol      Optional tolerance (mapped to `factr`).
  !! @param[in]    max_iter Optional maximum number of iterations.
  !! @param[in]    verbose  Optional verbosity flag.
  !!
  !! @details
  !! This uses the classic L-BFGS-B driver `setulb`. The objective function is
  !! computed by `compute_loss_and_gradient`.
  !!
  !! Training loss:
  !! $$ L = \sum_i \frac{1}{2} ||y_i - \hat{y}_i||^2 + \frac{\lambda}{2}||W||^2 $$
  !!
  !! @warning
  !! - Requires `training_inputs` and `training_targets` to be allocated.
  !! - Only one network may be trained at a time because of module-global state.
  subroutine nn_train_lbfgsb(this, tol, max_iter, verbose)
    class(type_nn), target, intent(inout) :: this
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    logical, intent(in), optional :: verbose

    integer :: n, m_lbfgs, iprint
    real(dp), allocatable :: x(:), l(:), u(:), g(:), wa(:)
    integer, allocatable :: nbd(:), iwa(:)
    real(dp) :: f, factr, pgtol
    real(dp), allocatable :: grad_work(:)
    character(len=60) :: task, csave
    logical :: lsave(4)
    integer :: isave(44)
    real(dp) :: dsave(29)
    integer :: max_it, iter
    logical :: verb

    external :: setulb

    if (.not. allocated(this%training_inputs)) then
      error stop "Error: Training data not set"
    end if
    if (callback_active) then
      error stop "Error: Cannot train multiple networks simultaneously"
    end if

    verb = .false.
    if (present(verbose)) verb = verbose

    factr = 1.0e7_dp
    if (present(tol)) factr = tol / epsilon(1.0_dp)

    pgtol = 1.0e-5_dp

    max_it = 1000
    if (present(max_iter)) max_it = max_iter

    n = this%n_params
    m_lbfgs = 10

    allocate(x(n), l(n), u(n), g(n), nbd(n))
    allocate(iwa(3*n))
    allocate(wa(2*m_lbfgs*n + 5*n + 11*m_lbfgs*m_lbfgs + 8*m_lbfgs))
    allocate(grad_work(n))  ! currently unused; kept for compatibility

    nbd = 0
    l = 0.0_dp
    u = 0.0_dp

    call this%network_to_vector(x)

    current_nn => this
    callback_active = .true.

    if (verb) then
      iprint = 1
      print *, "Starting L-BFGS-B optimization..."
      print *, "  Parameters (n):", n
      print *, "  Training samples:", this%n_samples
      print *, "  Max iterations:", max_it
    else
      iprint = -1
    end if

    task = 'START'
    iter = 0

    do while (task(1:2) == 'FG' .or. task == 'START' .or. task == 'NEW_X')
      call setulb(n, m_lbfgs, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)

      if (task(1:2) == 'FG') then
        call compute_loss_and_gradient(this, x, f, g)
      else if (task(1:5) == 'NEW_X') then
        iter = iter + 1
        if (iter > max_it) task = 'STOP: MAX ITERATIONS'
      else
        exit
      end if
    end do

    call this%vector_to_network(x)

    if (verb) then
      print *, ""
      print *, "=== L-BFGS-B Optimization Results ==="
      print *, "Status:", trim(task)
      print *, "Iterations:", iter
      print *, "Final loss:", f
      print *, "Final RMS error:", sqrt(2.0_dp * f / (this%n_samples * this%n_outputs))

      if (allocated(this%validation_inputs)) then
        ! NOTE: validation_inputs/targets may be normalized working copies,
        ! but evaluate() expects ORIGINAL space. In typical use, you'd pass
        ! validation_inputs_original/targets_original here.
        if (allocated(this%validation_inputs_original)) then
          print *, "Validation RMS:", this%evaluate(this%validation_inputs_original, this%validation_targets_original)
        else
          print *, "Validation RMS:", this%evaluate(this%validation_inputs, this%validation_targets)
        end if
      end if
    end if

    current_nn => null()
    callback_active = .false.
    deallocate(x, l, u, g, nbd, iwa, wa, grad_work)
  end subroutine nn_train_lbfgsb

  !> Compute loss and gradient for the current parameter vector (for L-BFGS-B).
  !!
  !! @param[inout] net    Network object (will be updated with `params`).
  !! @param[in]    params Parameter vector (length n_params).
  !! @param[out]   loss   Scalar objective value.
  !! @param[out]   grad   Gradient vector (length n_params).
  !!
  !! @details
  !! Performs full-batch backpropagation over all training samples.
  !!
  !! @note
  !! This uses `net%training_inputs` / `net%training_targets` (working copies),
  !! which may be normalized if `normalize_data()` has been called.
  subroutine compute_loss_and_gradient(net, params, loss, grad)
    type(type_nn), intent(inout) :: net
    real(dp), intent(in) :: params(:)
    real(dp), intent(out) :: loss, grad(:)

    integer :: i, j, k, idx
    real(dp) :: error, delta_h
    real(dp), allocatable :: delta_hidden(:), delta_output(:)

    call net%vector_to_network(params)

    loss = 0.0_dp
    grad = 0.0_dp
    allocate(delta_hidden(net%n_hidden))
    allocate(delta_output(net%n_outputs))

    do i = 1, net%n_samples
      call net%forward(net%training_inputs(:,i))

      do j = 1, net%n_outputs
        error = net%output(j) - net%training_targets(j,i)
        loss = loss + 0.5_dp * error * error
      end do

      do j = 1, net%n_outputs
        error = net%output(j) - net%training_targets(j,i)
        delta_output(j) = error * activation_derivative(net%output(j), net%output_activation)
      end do

      do k = 1, net%n_hidden
        delta_h = 0.0_dp
        do j = 1, net%n_outputs
          delta_h = delta_h + delta_output(j) * net%weights_ho(k,j)
        end do
        delta_hidden(k) = delta_h * activation_derivative(net%hidden(k), net%hidden_activation)
      end do

      idx = 1
      do k = 1, net%n_hidden
        do j = 1, net%n_inputs
          grad(idx) = grad(idx) + delta_hidden(k) * net%training_inputs(j,i)
          idx = idx + 1
        end do
      end do

      do j = 1, net%n_outputs
        do k = 1, net%n_hidden
          grad(idx) = grad(idx) + delta_output(j) * net%hidden(k)
          idx = idx + 1
        end do
      end do

      do k = 1, net%n_hidden
        grad(idx) = grad(idx) + delta_hidden(k)
        idx = idx + 1
      end do

      do j = 1, net%n_outputs
        grad(idx) = grad(idx) + delta_output(j)
        idx = idx + 1
      end do
    end do

    if (net%l2_lambda > 0.0_dp) then
      idx = net%n_inputs * net%n_hidden + net%n_hidden * net%n_outputs
      loss = loss + 0.5_dp * net%l2_lambda * sum(params(1:idx)**2)
      grad(1:idx) = grad(1:idx) + net%l2_lambda * params(1:idx)
    end if

    deallocate(delta_hidden, delta_output)
  end subroutine compute_loss_and_gradient

  !> Derivative of activation function with respect to its input.
  !!
  !! @param[in] y               Activation output value (not the pre-activation).
  !! @param[in] activation_type Activation ID.
  !! @return dy                 Derivative d(activation)/dz evaluated at the corresponding z.
  !!
  !! @details
  !! Uses output-based derivatives where available:
  !! - tanh: 1 - y^2
  !! - sigmoid: y(1-y)
  !!
  !! @warning
  !! For ReLU, the derivative depends on pre-activation z, but we approximate
  !! using the output y:
  !! - if y>0 => derivative 1
  !! - else => derivative 0
  function activation_derivative(y, activation_type) result(dy)
    real(dp), intent(in) :: y
    integer, intent(in) :: activation_type
    real(dp) :: dy

    select case(activation_type)
      case(ACTIVATION_TANH)
        dy = 1.0_dp - y * y
      case(ACTIVATION_SIGMOID)
        dy = y * (1.0_dp - y)
      case(ACTIVATION_RELU)
        if (y > 0.0_dp) then
          dy = 1.0_dp
        else
          dy = 0.0_dp
        end if
      case(ACTIVATION_LINEAR)
        dy = 1.0_dp
      case default
        dy = 1.0_dp
    end select
  end function activation_derivative

  !> Finite-difference gradient check (debugging utility).
  !!
  !! @param[inout] net     Network object.
  !! @param[in]    verbose Optional verbosity flag.
  !!
  !! @details
  !! Compares analytic gradients from backprop to numerical gradients computed
  !! via forward finite differences:
  !! $$g_i \approx (L(\theta_i+\epsilon)-L(\theta))/\epsilon$$
  !!
  !! @warning
  !! This is expensive: O(n_params) evaluations of the loss.
  subroutine nn_check_gradient(net, verbose)
    class(type_nn), intent(inout) :: net
    logical, intent(in), optional :: verbose

    real(dp), allocatable :: params(:), grad_analytic(:), grad_numeric(:), g_tmp(:)
    real(dp) :: loss1, loss2, eps, max_error, rel_error
    integer :: i, n_check
    logical :: verb

    verb = .false.
    if (present(verbose)) verb = verbose

    eps = 1.0e-5_dp
    n_check = min(10, net%n_params)

    allocate(params(net%n_params))
    allocate(grad_analytic(net%n_params))
    allocate(grad_numeric(net%n_params))
    allocate(g_tmp(net%n_params))

    call net%network_to_vector(params)

    call compute_loss_and_gradient(net, params, loss1, grad_analytic)

    do i = 1, net%n_params
      params(i) = params(i) + eps
      call compute_loss_and_gradient(net, params, loss2, g_tmp)
      params(i) = params(i) - eps
      grad_numeric(i) = (loss2 - loss1) / eps
    end do

    max_error = maxval(abs(grad_analytic - grad_numeric))
    rel_error = max_error / (maxval(abs(grad_analytic)) + 1.0e-8_dp)

    if (verb) then
      print *, ""
      print *, "=== Gradient Check ==="
      print *, "Max absolute error:", max_error
      print *, "Max relative error:", rel_error
      print *, ""
      print *, "Sample gradients (first", n_check, "parameters):"
      do i = 1, n_check
        print '(A,I3,A,E12.4,A,E12.4,A,E12.4)', &
          "  Param", i, ": Analytic=", grad_analytic(i), &
          " Numeric=", grad_numeric(i), &
          " Error=", abs(grad_analytic(i) - grad_numeric(i))
      end do
      print *, ""
      if (rel_error < 1.0e-5_dp) then
        print *, "Gradient check PASSED (error < 1e-5)"
      else if (rel_error < 1.0e-3_dp) then
        print *, "Gradient check OK (error < 1e-3, but could be better)"
      else
        print *, "Gradient check FAILED (error too large!)"
      end if
      print *, ""
    end if

    deallocate(params, grad_analytic, grad_numeric, g_tmp)
  end subroutine nn_check_gradient

end module module_nn
